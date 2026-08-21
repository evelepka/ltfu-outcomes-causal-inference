# =============================================================================
#  CLONE-CENSOR-WEIGHT (CCW) — ANÁLISE DO EFEITO DA LTFU NA MORTALIDADE
#  Versão comentada, passo a passo, para aprendizado.
#
#  Autoria da análise: Evelyn Lepka de Lima (com apoio de codificação)
#
#  ---------------------------------------------------------------------------
#  A PERGUNTA CAUSAL
#  ---------------------------------------------------------------------------
#  "Qual o efeito de desvincular-se do cuidado (LTFU) durante o tratamento de
#   tuberculose, em comparação com permanecer engajado, sobre a mortalidade?"
#
#  O DESAFIO:
#  A LTFU é uma "estratégia sustentada" que só se revela DEPOIS do início do
#  tratamento. No tempo zero (início do tratamento) não sabemos quem vai
#  desvincular. Além disso, para ser CLASSIFICADO como LTFU a pessoa precisa
#  SOBREVIVER até desvincular (a regra dos 30 dias). Isso cria dois problemas:
#
#   (1) IMMORTAL TIME BIAS: quem morre cedo, ainda em tratamento, nunca entra
#       no grupo LTFU. Logo o grupo LTFU parece "protegido" no início — um
#       artefato, não um efeito real.
#   (2) CONFUNDIMENTO: quem desvincula difere de quem fica (idade, vulnerabilidade).
#
#  COMO O CCW RESOLVE:
#   - CLONAGEM no tempo zero: cada pessoa entra nos DOIS braços. Assim, mortes
#     precoces (antes das estratégias divergirem) contam nos dois grupos —
#     como numa análise por intenção-de-tratar (ITT) de um ensaio clínico.
#   - CENSURA ao desviar: um clone é censurado quando o comportamento real da
#     pessoa desvia da estratégia atribuída.
#   - IPCW (inverse-probability-of-censoring weights): repondera para corrigir
#     a seleção causada por essa censura artificial.
# =============================================================================

import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import OneHotEncoder
from sklearn.compose import ColumnTransformer
import warnings
warnings.filterwarnings("ignore")

# Constantes do estudo
MONTH   = 30.4     # dias por mês (média) — para converter dias em meses
HORIZON = 24       # horizonte de seguimento em meses (nosso desfecho é 24 meses)

# Covariáveis basais usadas para AJUSTE (no modelo de pesos IPCW).
# São medidas no início do tratamento, então não são afetadas pela exposição.
COVS = ['age_group','sex','hiv_aids','homelessness','hosp_admission',
        'drug_use','clinical_clean','dot_status']


# =============================================================================
#  PASSO 1 — CONSTRUIR A LINHA DO TEMPO DE CADA PESSOA
#  ---------------------------------------------------------------------------
#  Precisamos, para cada indivíduo, de três "relógios" medidos a partir do
#  INÍCIO DO TRATAMENTO (o nosso tempo zero):
#     - quando desvinculou (se LTFU)   -> m_dis
#     - quando o tratamento terminou   -> m_end   (conclusão ou fechamento)
#     - quando morreu (se morreu)      -> m_death
#
#  NOTA IMPORTANTE sobre a data de LTFU:
#  A regra brasileira classifica LTFU após 30 dias sem tratamento. A data de
#  fechamento do caso (end_date) é ~30 dias APÓS a interrupção. Então a data
#  real de interrupção = end_date - 30 dias. É isso que usamos como m_dis.
# =============================================================================

def preparar_timeline(caminho_cohort):
    df = pd.read_csv(caminho_cohort, low_memory=False)

    # Converter as três datas em datetime
    for c in ['best_start', 'end_date', 'death_date']:
        df[c] = pd.to_datetime(df[c], errors='coerce')

    # Marcar quem é LTFU
    df['is_ltfu'] = (df['itt_group'] == 'Loss to follow-up')

    # Data de desvinculação = end_date - 30 dias (só para LTFU)
    df['ltfu_date'] = pd.NaT
    df.loc[df['is_ltfu'], 'ltfu_date'] = df.loc[df['is_ltfu'], 'end_date'] - pd.Timedelta(days=30)

    # Três relógios em DIAS a partir do início do tratamento
    df['t_disengage'] = (df['ltfu_date']  - df['best_start']).dt.days
    df['t_txend']     = (df['end_date']   - df['best_start']).dt.days
    df['t_death']     = (df['death_date'] - df['best_start']).dt.days
    df['died']        = (df['event_d'] == 1).astype(int)

    # Converter os relógios de DIAS para MESES (inteiros), limitando ao horizonte.
    # np.maximum(...,0) evita valores negativos (raros casos com datas invertidas).
    df['m_dis']   = np.floor(np.maximum(df['t_disengage'], 0) / MONTH).clip(upper=HORIZON)
    df['m_end']   = np.floor(np.maximum(df['t_txend'],     0) / MONTH).clip(upper=HORIZON)
    df['m_death'] = np.where(df['died'] == 1,
                             np.floor(np.maximum(df['t_death'], 0) / MONTH),
                             np.nan)   # NaN quando não houve morte
    return df


# =============================================================================
#  PASSO 2 — CLONAR E CENSURAR (o coração do método)
#  ---------------------------------------------------------------------------
#  Criamos um "braço" (arm) atribuindo TODOS os indivíduos a uma estratégia.
#  A função abaixo constrói UM braço:
#     kind = 'remain'          -> estratégia "permanecer engajado"
#     kind = ('janela', meses) -> estratégia "desvincular na janela de meses"
#
#  A LÓGICA DA CENSURA (deviation):
#  Um clone é censurado no momento em que a pessoa DESVIA da estratégia:
#   - No braço "permanecer": desvia quando desvincula (LTFU) -> censura em m_dis
#   - No braço "desvincular": desvia quando (a) completa o tratamento engajado
#     [censura em m_end] OU (b) desvincula fora da janela alvo [censura em m_dis]
#
#  A REGRA DE OURO (por que mortes precoces contam nos dois braços):
#  A pessoa só "desvia" se estiver VIVA no momento do desvio. Se ela MORRE antes
#  de desviar, a morte é contada (não há censura). Como no tempo zero todos são
#  compatíveis com as duas estratégias, uma morte precoce (ainda em tratamento)
#  aparece nos DOIS braços. É exatamente isso que corrige o immortal time bias.
# =============================================================================

def construir_braco(dat, kind):
    d = dat.copy()
    md = d['m_death']
    # helper: a pessoa está VIVA no mês do evento de desvio?
    vivo_em = lambda mes_evento: md.isna() | (md > mes_evento)

    if kind == 'remain':
        # Braço "permanecer engajado":
        #   - LTFU desviam ao desvincular (censura em m_dis), SE vivos lá
        #   - não-LTFU aderem sempre (nunca desviam) -> m_dev = infinito
        d['m_dev'] = np.where(d['is_ltfu'],
                              np.where(vivo_em(d['m_dis']), d['m_dis'], np.inf),
                              np.inf)
    else:
        # Braço "desvincular na janela alvo":
        meses_alvo = kind[1]
        desvinc_na_janela = d['is_ltfu'] & d['m_dis'].isin(meses_alvo)
        # desvio por completar o tratamento engajado (se vivo em m_end)
        dev_completa = np.where(vivo_em(d['m_end']), d['m_end'], np.inf)
        # desvio por desvincular FORA da janela (se vivo em m_dis)
        dev_janela_errada = np.where(d['is_ltfu'] & ~desvinc_na_janela & vivo_em(d['m_dis']),
                                     d['m_dis'], np.inf)
        m_dev = np.minimum(dev_completa, dev_janela_errada)
        # quem ADERE à estratégia (desvinculou na janela certa) nunca é censurado
        d['m_dev'] = np.where(desvinc_na_janela, np.inf, m_dev)

    # Fim do seguimento deste clone = min(mês do desvio, horizonte)
    d['stop'] = np.minimum(d['m_dev'], HORIZON)
    # Mês do evento (morte) SÓ conta se ocorrer ANTES do desvio/censura
    d['event_month'] = np.where(md < d['stop'], md, np.inf)
    # Número de meses que o clone contribui (inclui o mês da morte, se houver)
    nm = np.where(d['event_month'] < np.inf, d['event_month'] + 1, d['stop'])
    d['nmonths'] = np.clip(nm, 0, None).astype(int)
    return d[d['nmonths'] > 0]


# =============================================================================
#  PASSO 3 — EXPANDIR PARA "PESSOA-MÊS" (formato longo)
#  ---------------------------------------------------------------------------
#  Modelos de tempo discreto trabalham com uma linha por pessoa por mês.
#  Cada clone que segue N meses vira N linhas (mês 0, 1, 2, ...).
#  Marcamos em cada linha:
#     death      = 1 se a morte ocorre NAQUELE mês
#     dev_next   = 1 se o clone será CENSURADO (desvio) ao final daquele mês
#                  (isso alimenta o modelo de pesos IPCW)
# =============================================================================

def expandir(clones):
    n = clones['nmonths'].values.astype(int)
    idx = np.repeat(np.arange(len(clones)), n)          # repete cada linha N vezes
    L = clones.iloc[idx].reset_index(drop=True)
    L['month'] = np.concatenate([np.arange(k) for k in n])   # 0,1,2,... por clone
    L['death'] = ((L['event_month'] < np.inf) & (L['month'] == L['event_month'])).astype(int)
    L['dev_next'] = ((L['m_dev'] < HORIZON) &
                     (L['month'] == L['stop'] - 1) &
                     (L['death'] == 0)).astype(int)
    L['rid'] = idx    # identificador do clone (para o produto cumulativo dos pesos)
    return L


# =============================================================================
#  PASSO 4 — CALCULAR OS PESOS IPCW (a correção da censura)
#  ---------------------------------------------------------------------------
#  A censura artificial (no desvio) é INFORMATIVA: quem desvia difere de quem
#  adere. Se ignorarmos, os grupos ficam enviesados. O IPCW corrige repondendo
#  cada pessoa-mês pelo INVERSO da probabilidade de "sobreviver à censura".
#
#  Como se calcula:
#   1) Modelo logístico da probabilidade de DESVIO em cada mês, dado covariáveis
#      basais + mês (denominador).
#   2) Modelo só com o mês (numerador) -> gera pesos ESTABILIZADOS (mais estáveis).
#   3) Peso = produto cumulativo, ao longo dos meses, de (1 - p_num)/(1 - p_den).
#
#  Pesos estabilizados devem ter média ~1. Truncamos extremos (0,5% e 99,5%)
#  para evitar que poucos indivíduos dominem a estimativa.
# =============================================================================

def calcular_ipcw(L):
    for c in COVS:
        L[c] = L[c].astype(str)
    L['mS'] = L['month'].astype(str)

    if L['dev_next'].sum() > 0:
        # codificar covariáveis + mês como variáveis dummy (one-hot)
        enc  = ColumnTransformer([('oh', OneHotEncoder(handle_unknown='ignore',
                                    sparse_output=True), COVS + ['mS'])])
        encm = ColumnTransformer([('oh', OneHotEncoder(handle_unknown='ignore',
                                    sparse_output=True), ['mS'])])
        Xden = enc.fit_transform(L)     # matriz do denominador (com covariáveis)
        Xnum = encm.fit_transform(L)    # matriz do numerador (só mês)
        y = L['dev_next'].values

        # probabilidade de desvio (denominador e numerador)
        p_den = LogisticRegression(max_iter=150, C=1e6).fit(Xden, y).predict_proba(Xden)[:, 1]
        p_num = LogisticRegression(max_iter=150, C=1e6).fit(Xnum, y).predict_proba(Xnum)[:, 1]

        # probabilidade de NÃO desviar naquele mês
        t = pd.DataFrame({'g': L['rid'].values, 'd': 1 - p_den, 'n': 1 - p_num})
        # produto cumulativo ao longo dos meses de cada clone
        t['cum_den'] = t.groupby('g')['d'].cumprod()
        t['cum_num'] = t.groupby('g')['n'].cumprod()
        L['sw'] = (t['cum_num'] / t['cum_den']).values   # peso estabilizado
    else:
        L['sw'] = 1.0

    # truncar pesos extremos
    lo, hi = L['sw'].quantile([0.005, 0.995])
    L['swt'] = L['sw'].clip(lo, hi)
    return L


# =============================================================================
#  PASSO 5 — INCIDÊNCIA CUMULATIVA PONDERADA (o desfecho)
#  ---------------------------------------------------------------------------
#  Com os pesos, estimamos o risco de morte acumulado ao longo dos meses:
#     - risco (hazard) de cada mês = média PONDERADA das mortes daquele mês
#     - sobrevivência = produto de (1 - hazard) ao longo dos meses
#     - incidência cumulativa = 1 - sobrevivência
# =============================================================================

def risco_ponderado(L):
    hazard = L.groupby('month').apply(lambda g: np.average(g['death'], weights=g['swt']))
    incidencia = 1 - (1 - hazard).cumprod()
    return incidencia


# =============================================================================
#  PASSO 6 — UM CONTRASTE COMPLETO (desvincular vs permanecer)
#  ---------------------------------------------------------------------------
#  Junta tudo: constrói os dois braços, calcula IPCW, estima o risco de cada um
#  e devolve as medidas de efeito (razão de risco RR e diferença de risco RD)
#  em 6 e 24 meses.
# =============================================================================

def contraste(dat, meses_janela, rotulo=""):
    # braço "desvincular" (arm=1) e braço "permanecer" (arm=0)
    dis = expandir(construir_braco(dat, ('janela', meses_janela))); dis['cid'] = 'D_' + dis['rid'].astype(str)
    rem = expandir(construir_braco(dat, 'remain'));                 rem['cid'] = 'R_' + rem['rid'].astype(str)
    L = pd.concat([dis.assign(A=1), rem.assign(A=0)], ignore_index=True)

    # IPCW por braço (usa 'cid' como identificador único do clone)
    L['rid'] = L['cid']  # reaproveita a coluna de agrupamento
    L = calcular_ipcw(L)

    ci1 = risco_ponderado(L[L.A == 1])   # risco no braço desvincular
    ci0 = risco_ponderado(L[L.A == 0])   # risco no braço permanecer

    def em(ci, mes):
        s = ci.loc[ci.index <= mes]
        return s.iloc[-1] if len(s) else np.nan

    r6_1, r6_0   = em(ci1, 5), em(ci0, 5)     # risco em ~6 meses (meses 0..5)
    r24_1, r24_0 = ci1.iloc[-1], ci0.iloc[-1] # risco em 24 meses

    print(f"### {rotulo} ###")
    print(f"  6 meses : desvincular {r6_1*100:5.2f}%  vs  permanecer {r6_0*100:5.2f}%"
          f"   RR={r6_1/r6_0:.2f}  RD={100*(r6_1-r6_0):+.2f}pp")
    print(f"  24 meses: desvincular {r24_1*100:5.2f}%  vs  permanecer {r24_0*100:5.2f}%"
          f"   RR={r24_1/r24_0:.2f}  RD={100*(r24_1-r24_0):+.2f}pp\n")
    return dict(rr24=r24_1/r24_0, rd24=100*(r24_1-r24_0), rr6=r6_1/r6_0)


# =============================================================================
#  EXECUÇÃO PRINCIPAL
# =============================================================================
if __name__ == "__main__":
    ROOT = Path.home() / ("Library/CloudStorage/GoogleDrive-evelynlepka@gmail.com/"
                          "My Drive/Abandonment Outcomes/Abandonment Paper")
    tl = preparar_timeline(ROOT / "ITT_Analysis/data/itt_cohort.csv")

    # As três estratégias de timing (sua "opção 2": três níveis)
    print("="*64)
    print("EFEITO DA LTFU NA MORTALIDADE — por timing da desvinculação")
    print("="*64 + "\n")
    contraste(tl, [0],          "Abandono primário (mês 0)")
    contraste(tl, [1, 2],       "Fase intensiva (meses 1-2)")
    contraste(tl, [3, 4, 5, 6], "Fase de continuação (meses 3-6)")

    # Gradiente mês a mês (figura secundária/exploratória)
    print("="*64)
    print("GRADIENTE mês a mês (exploratório — cauda tardia é instável)")
    print("="*64 + "\n")
    for m in range(1, 7):
        contraste(tl, [m], f"Desvinculação no mês {m}")

# =============================================================================
#  NOTA SOBRE OS INTERVALOS DE CONFIANÇA
#  ---------------------------------------------------------------------------
#  Os ICs vêm de BOOTSTRAP: reamostramos indivíduos com reposição centenas de
#  vezes, refazemos TODO o pipeline (clonar, censurar, IPCW, estimar) em cada
#  reamostra, e tomamos os percentis 2,5 e 97,5 das estimativas. Isso captura
#  corretamente a incerteza, inclusive a introduzida pela estimação dos pesos.
#  (O script de bootstrap é separado por ser computacionalmente pesado.)
# =============================================================================
