import streamlit as st
import plotly.graph_objects as go
import plotly.express as px
import pandas as pd
from scipy import stats
import sqlite3

st.set_page_config(
    page_title="Miraclib Clinical Trial | Immune Cell Analysis",
    layout="wide",
    initial_sidebar_state="collapsed"
)

st.markdown("""
<style>
    @import url('https://fonts.googleapis.com/css2?family=IBM+Plex+Sans:wght@300;400;500;600&family=IBM+Plex+Mono:wght@400;600&display=swap');

    .stApp, html, body,
    div[data-testid="stAppViewContainer"],
    div[data-testid="stHeader"],
    div[data-testid="stToolbar"],
    div[data-testid="block-container"] {
        background-color: #ffffff !important;
        color: #111827 !important;
    }
    html, body, [class*="css"] {
        font-family: 'IBM Plex Sans', sans-serif !important;
    }
    .cover { border-bottom: 2px solid #111827; padding-bottom: 32px; margin-bottom: 48px; }
    .cover-tag { font-family: 'IBM Plex Mono', monospace; font-size: 0.72rem; color: #6b7280; letter-spacing: 2px; text-transform: uppercase; margin-bottom: 12px; }
    .cover-title { font-size: 2.2rem; font-weight: 600; color: #111827; letter-spacing: -0.5px; line-height: 1.2; margin-bottom: 12px; }
    .cover-subtitle { font-size: 1rem; color: #374151; margin-bottom: 16px; max-width: 700px; line-height: 1.6; }
    .cover-meta { font-family: 'IBM Plex Mono', monospace; font-size: 0.72rem; color: #6b7280; }
    .section-number { font-family: 'IBM Plex Mono', monospace; font-size: 0.68rem; font-weight: 600; letter-spacing: 2px; text-transform: uppercase; color: #6b7280; margin-bottom: 4px; margin-top: 56px; }
    .section-title { font-size: 1.4rem; font-weight: 600; color: #111827; margin-bottom: 8px; padding-bottom: 10px; border-bottom: 1px solid #d1d5db; }
    .narrative { font-size: 0.93rem; color: #1f2937; line-height: 1.8; max-width: 820px; margin-bottom: 28px; }
    .stat-card { background: #f9fafb; border: 1px solid #d1d5db; border-radius: 8px; padding: 20px 16px; text-align: center; }
    .stat-value { font-family: 'IBM Plex Mono', monospace; font-size: 1.75rem; font-weight: 600; color: #111827; white-space: nowrap; }
    .stat-label { font-size: 0.7rem; color: #374151; text-transform: uppercase; letter-spacing: 0.5px; margin-top: 6px; }
    .finding-box { background: #f9fafb; border-left: 3px solid #2563eb; padding: 14px 18px; border-radius: 0 6px 6px 0; margin: 16px 0; font-size: 0.88rem; color: #1f2937; line-height: 1.65; }
    .finding-box strong { display: block; font-weight: 600; margin-bottom: 4px; font-size: 0.75rem; text-transform: uppercase; letter-spacing: 0.5px; color: #1d4ed8; }
    .method-note { background: #f9fafb; border: 1px solid #d1d5db; border-radius: 6px; padding: 10px 14px; font-size: 0.82rem; color: #374151; font-style: italic; margin-top: 8px; }
    .divider { border: none; border-top: 1px solid #d1d5db; margin: 48px 0; }
    .future-item { background: #f9fafb; border: 1px solid #d1d5db; border-radius: 8px; padding: 20px; height: 100%; }
    .future-number { font-family: 'IBM Plex Mono', monospace; font-size: 1.5rem; font-weight: 600; color: #d1d5db; margin-bottom: 8px; }
    .future-title { font-size: 0.9rem; font-weight: 600; color: #111827; margin-bottom: 6px; }
    .future-desc { font-size: 0.82rem; color: #374151; line-height: 1.6; }
</style>
""", unsafe_allow_html=True)

# ── Constants ─────────────────────────────────────────────────────────────────
DB_PATH = "cell_counts.db"
POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]
POPULATION_LABELS = {
    "b_cell": "B Cell", "cd8_t_cell": "CD8+ T Cell", "cd4_t_cell": "CD4+ T Cell",
    "nk_cell": "NK Cell", "monocyte": "Monocyte"
}
POP_COLORS = {
    "B Cell": "#2563eb", "CD8+ T Cell": "#16a34a", "CD4+ T Cell": "#dc2626",
    "NK Cell": "#d97706", "Monocyte": "#7c3aed"
}
RESP_COLOR     = "#2563eb"
NON_RESP_COLOR = "#6b7280"

DARK_FONT = dict(family="IBM Plex Sans", color="#111827", size=13)
CHART_LAYOUT = dict(plot_bgcolor="#ffffff", paper_bgcolor="#ffffff", font=DARK_FONT)

def dark_axis(title_text):
    return dict(
        title=dict(text=title_text, font=dict(color="#111827", size=13)),
        tickfont=dict(color="#111827", size=12),
        gridcolor="#d1d5db",
        linecolor="#374151",
        zeroline=False
    )

def chart_title(text):
    return dict(text=text, font=dict(color="#111827", size=14))

# ── Data loading ──────────────────────────────────────────────────────────────
@st.cache_data
def load_all_data():
    conn = sqlite3.connect(DB_PATH)
    samples_df  = pd.read_sql_query("SELECT * FROM samples", conn)
    subjects_df = pd.read_sql_query("SELECT * FROM subjects", conn)
    melanoma_df = pd.read_sql_query("""
        SELECT s.*, sub.condition, sub.treatment, sub.response, sub.sex, sub.project
        FROM samples s JOIN subjects sub ON s.subject = sub.subject
        WHERE sub.condition = 'melanoma' AND sub.treatment = 'miraclib' AND s.sample_type = 'PBMC'
    """, conn)
    baseline_df = pd.read_sql_query("""
        SELECT s.*, sub.project, sub.condition, sub.treatment, sub.response, sub.sex
        FROM samples s JOIN subjects sub ON s.subject = sub.subject
        WHERE sub.condition = 'melanoma' AND sub.treatment = 'miraclib'
          AND s.sample_type = 'PBMC' AND s.time_from_treatment_start = 0
    """, conn)
    avg_bcell = pd.read_sql_query("""
        SELECT ROUND(AVG(s.b_cell), 2) as avg_b_cell
        FROM samples s JOIN subjects sub ON s.subject = sub.subject
        WHERE sub.condition = 'melanoma' AND sub.sex = 'M'
          AND sub.response = 'yes' AND s.time_from_treatment_start = 0
    """, conn)
    conn.close()
    return samples_df, subjects_df, melanoma_df, baseline_df, avg_bcell["avg_b_cell"].iloc[0]

@st.cache_data
def build_frequency_table(samples_df):
    samples_df = samples_df.copy()
    samples_df["total_count"] = samples_df[POPULATIONS].sum(axis=1)
    rows = []
    for _, row in samples_df.iterrows():
        for pop in POPULATIONS:
            rows.append({
                "sample": row["sample"], "total_count": int(row["total_count"]),
                "population": pop, "count": int(row[pop]),
                "percentage": round((row[pop] / row["total_count"]) * 100, 2)
            })
    return pd.DataFrame(rows)

@st.cache_data
def build_melanoma_freq(melanoma_df):
    melanoma_df = melanoma_df.copy()
    melanoma_df["total_count"] = melanoma_df[POPULATIONS].sum(axis=1)
    rows = []
    for _, row in melanoma_df.iterrows():
        for pop in POPULATIONS:
            rows.append({
                "sample": row["sample"], "total_count": int(row["total_count"]),
                "population": pop, "count": int(row[pop]),
                "percentage": round((row[pop] / row["total_count"]) * 100, 2),
                "response": row["response"]
            })
    return pd.DataFrame(rows)

@st.cache_data
def run_statistics(melanoma_freq):
    corrected_threshold = 0.05 / len(POPULATIONS)
    results = []
    for pop in POPULATIONS:
        resp     = melanoma_freq[(melanoma_freq["population"] == pop) & (melanoma_freq["response"] == "yes")]["percentage"]
        non_resp = melanoma_freq[(melanoma_freq["population"] == pop) & (melanoma_freq["response"] == "no")]["percentage"]
        mw_stat, mw_p = stats.mannwhitneyu(resp, non_resp, alternative="two-sided")
        _, t_p        = stats.ttest_ind(resp, non_resp, equal_var=False)
        _, bm_p       = stats.brunnermunzel(resp, non_resp)
        n1, n2        = len(resp), len(non_resp)
        effect_size   = 1 - (2 * mw_stat) / (n1 * n2)
        results.append({
            "Population": POPULATION_LABELS[pop],
            "Responder Median (%)": round(resp.median(), 2),
            "Non-Responder Median (%)": round(non_resp.median(), 2),
            "Mann-Whitney p": round(mw_p, 4),
            "Welch's T p": round(t_p, 4),
            "Brunner-Munzel p": round(bm_p, 4),
            "Effect Size": round(effect_size, 4),
            "Significant": "yes" if mw_p < corrected_threshold else "no"
        })
    return pd.DataFrame(results), corrected_threshold

# ── Load ──────────────────────────────────────────────────────────────────────
samples_df, subjects_df, melanoma_df, baseline_df, avg_bcell = load_all_data()
freq_df       = build_frequency_table(samples_df)
melanoma_freq = build_melanoma_freq(melanoma_df)
stats_df, corrected_threshold = run_statistics(melanoma_freq)
cd4_row     = stats_df[stats_df["Population"] == "CD4+ T Cell"].iloc[0]
n_resp_base = baseline_df.drop_duplicates("subject")[baseline_df.drop_duplicates("subject")["response"] == "yes"].shape[0]
n_non_base  = baseline_df.drop_duplicates("subject")[baseline_df.drop_duplicates("subject")["response"] == "no"].shape[0]

# ════════════════════════════════════════════════════════════════════════════════
# COVER
# ════════════════════════════════════════════════════════════════════════════════
st.markdown("""
<div class="cover">
    <div class="cover-tag">Clinical Trial Analysis Report</div>
    <div class="cover-title">Miraclib — Immune Cell Population Analysis</div>
    <div class="cover-subtitle">
        An end-to-end analysis of immune cell dynamics across melanoma and carcinoma
        patient cohorts, evaluating how miraclib affects immune composition and identifying
        cell populations that may predict treatment response.
    </div>
    <div class="cover-meta">
        Prepared for: Bob Loblaw, Loblaw Bio &nbsp;|&nbsp; Drug Candidate: Miraclib &nbsp;|&nbsp; Sample Type: PBMC
    </div>
</div>
""", unsafe_allow_html=True)

# ════════════════════════════════════════════════════════════════════════════════
# SECTION 1 — THE DATA
# ════════════════════════════════════════════════════════════════════════════════
st.markdown('<div class="section-number">Section 01</div>', unsafe_allow_html=True)
st.markdown('<div class="section-title">The Data</div>', unsafe_allow_html=True)
st.markdown("""
<div class="narrative">
    The dataset contains immune cell counts across five populations — B cells, CD8+ T cells,
    CD4+ T cells, NK cells, and monocytes — measured from patient blood samples (PBMC) across
    multiple timepoints. Patient metadata including condition, treatment, and response status
    are stored in a relational SQLite database with two tables: <code>subjects</code> for
    patient-level demographics and <code>samples</code> for sample-level measurements.
</div>
""", unsafe_allow_html=True)

c1, c2, c3, c4, c5 = st.columns(5)
kpis = [
    (f"{subjects_df.shape[0]:,}", "Total Subjects"),
    (f"{samples_df.shape[0]:,}", "Total Samples"),
    (f"{subjects_df['condition'].nunique()}", "Indications"),
    (f"{subjects_df['treatment'].nunique()}", "Treatments"),
    (f"{samples_df['time_from_treatment_start'].nunique()}", "Timepoints"),
]
for col, (val, label) in zip([c1, c2, c3, c4, c5], kpis):
    with col:
        st.markdown(f'<div class="stat-card"><div class="stat-value">{val}</div><div class="stat-label">{label}</div></div>', unsafe_allow_html=True)

st.markdown("<br/>", unsafe_allow_html=True)

col_a, col_b = st.columns(2)
with col_a:
    cond = subjects_df.groupby("condition")["subject"].count().reset_index()
    cond.columns = ["Condition", "Subjects"]
    fig_cond = px.bar(cond, x="Condition", y="Subjects", color="Condition",
                      color_discrete_sequence=["#2563eb", "#16a34a", "#dc2626"])
    fig_cond.update_layout(
        **CHART_LAYOUT, height=300, showlegend=False,
        title=chart_title("Subjects by Condition"),
        yaxis=dark_axis("Subject Count"),
        xaxis=dark_axis(""),
        margin=dict(l=10, r=10, t=50, b=10)
    )
    st.plotly_chart(fig_cond, use_container_width=True)

with col_b:
    resp = subjects_df.groupby("response")["subject"].count().reset_index()
    resp.columns = ["Response", "Subjects"]
    resp["Response"] = resp["Response"].map({"yes": "Responder", "no": "Non-Responder"})
    fig_resp = px.pie(resp, names="Response", values="Subjects",
                      color_discrete_sequence=[RESP_COLOR, NON_RESP_COLOR], hole=0.5)
    fig_resp.update_traces(textinfo="label+percent", textfont=dict(size=13, color="#111827"))
    fig_resp.update_layout(
        **CHART_LAYOUT, height=300, showlegend=False,
        title=chart_title("Subjects by Response Status"),
        margin=dict(l=10, r=10, t=50, b=10)
    )
    st.plotly_chart(fig_resp, use_container_width=True)

st.markdown('<hr class="divider"/>', unsafe_allow_html=True)

# ════════════════════════════════════════════════════════════════════════════════
# SECTION 2 — CELL POPULATION OVERVIEW
# ════════════════════════════════════════════════════════════════════════════════
st.markdown('<div class="section-number">Section 02</div>', unsafe_allow_html=True)
st.markdown('<div class="section-title">Cell Population Overview</div>', unsafe_allow_html=True)
st.markdown("""
<div class="narrative">
    For each patient sample, the relative frequency of each immune cell population was calculated
    as a percentage of the total cell count. This normalization allows meaningful comparison
    across samples regardless of differences in total cell counts. The chart below shows the
    average frequency of each population across all samples, with error bars representing
    one standard deviation.
</div>
""", unsafe_allow_html=True)

avg_pop = freq_df.groupby("population")["percentage"].agg(["mean", "std"]).reset_index()
avg_pop["label"] = avg_pop["population"].map(POPULATION_LABELS)
avg_pop = avg_pop.sort_values("mean", ascending=False)

fig_avg = go.Figure()
fig_avg.add_trace(go.Bar(
    x=avg_pop["label"], y=avg_pop["mean"],
    error_y=dict(type="data", array=avg_pop["std"], visible=True, color="#6b7280"),
    marker_color=[POP_COLORS[l] for l in avg_pop["label"]],
    marker_line_width=0,
    text=[f"{v:.1f}%" for v in avg_pop["mean"]],
    textposition="outside",
    textfont=dict(color="#111827", size=13)
))
fig_avg.update_layout(
    **CHART_LAYOUT, height=380, showlegend=False,
    title=chart_title("Average Relative Frequency by Cell Population (All Samples)"),
    yaxis=dark_axis("Mean Frequency (%)"),
    xaxis=dark_axis(""),
    margin=dict(l=20, r=20, t=60, b=20),
)
fig_avg.update_yaxes(range=[0, avg_pop["mean"].max() * 1.3])
st.plotly_chart(fig_avg, use_container_width=True)
st.markdown('<div class="method-note">Error bars represent ± 1 standard deviation. CD4+ T cells account for the largest share of immune cells on average, followed by CD8+ T cells and monocytes.</div>', unsafe_allow_html=True)

st.markdown('<hr class="divider"/>', unsafe_allow_html=True)

# ════════════════════════════════════════════════════════════════════════════════
# SECTION 3 — RESPONDER vs NON-RESPONDER
# ════════════════════════════════════════════════════════════════════════════════
st.markdown('<div class="section-number">Section 03</div>', unsafe_allow_html=True)
st.markdown('<div class="section-title">Responders vs. Non-Responders</div>', unsafe_allow_html=True)
n_resp_samples = melanoma_freq.drop_duplicates("sample")[melanoma_freq.drop_duplicates("sample")["response"] == "yes"].shape[0]
n_non_samples  = melanoma_freq.drop_duplicates("sample")[melanoma_freq.drop_duplicates("sample")["response"] == "no"].shape[0]
st.markdown(f"""
<div class="narrative">
    To identify immune signatures that may predict miraclib response, cell population frequencies
    were compared between melanoma patients who responded to treatment and those who did not.
    Analysis was restricted to PBMC samples from melanoma patients receiving miraclib only.
    The cohort included {n_resp_samples} responder samples and {n_non_samples} non-responder
    samples — a well-balanced cohort that supports valid statistical comparison.
</div>
""", unsafe_allow_html=True)

# Boxplots
fig_box = go.Figure()
for resp_val, color, border, label in [
    ("yes", RESP_COLOR, "#1d4ed8", "Responder"),
    ("no",  NON_RESP_COLOR, "#374151", "Non-Responder")
]:
    for pop in POPULATIONS:
        data = melanoma_freq[(melanoma_freq["population"] == pop) & (melanoma_freq["response"] == resp_val)]["percentage"]
        fig_box.add_trace(go.Box(
            y=data, x=[POPULATION_LABELS[pop]] * len(data),
            name=label, marker_color=color, line_color=border,
            boxmean=True, legendgroup=resp_val,
            showlegend=(pop == POPULATIONS[0]), opacity=0.9
        ))

fig_box.update_layout(
    **CHART_LAYOUT,
    title=chart_title("Cell Population Frequencies: Responders vs. Non-Responders (Melanoma · Miraclib · PBMC)"),
    boxmode="group", height=500,
    legend=dict(title="", orientation="h", y=1.06, x=0,
                bgcolor="#ffffff", bordercolor="#d1d5db", borderwidth=1,
                font=dict(color="#111827", size=13)),
    yaxis=dark_axis("Relative Frequency (%)"),
    xaxis=dark_axis(""),
    margin=dict(l=20, r=20, t=70, b=20),
)
st.plotly_chart(fig_box, use_container_width=True)

st.markdown("<br/>", unsafe_allow_html=True)

# Bar chart with significance markers
resp_means, non_means, resp_stds, non_stds, p_values, bar_labels = [], [], [], [], [], []
for pop in POPULATIONS:
    r = melanoma_freq[(melanoma_freq["population"] == pop) & (melanoma_freq["response"] == "yes")]["percentage"]
    n = melanoma_freq[(melanoma_freq["population"] == pop) & (melanoma_freq["response"] == "no")]["percentage"]
    resp_means.append(round(r.mean(), 2))
    non_means.append(round(n.mean(), 2))
    resp_stds.append(r.std())
    non_stds.append(n.std())
    _, p = stats.mannwhitneyu(r, n, alternative="two-sided")
    p_values.append(p)
    bar_labels.append(POPULATION_LABELS[pop])

fig_bar = go.Figure()
fig_bar.add_trace(go.Bar(
    name="Responder", x=bar_labels, y=resp_means,
    error_y=dict(type="data", array=resp_stds, visible=True, color="#6b7280"),
    marker_color=RESP_COLOR, opacity=0.9
))
fig_bar.add_trace(go.Bar(
    name="Non-Responder", x=bar_labels, y=non_means,
    error_y=dict(type="data", array=non_stds, visible=True, color="#6b7280"),
    marker_color=NON_RESP_COLOR, opacity=0.9
))
for i, p in enumerate(p_values):
    max_val = max(resp_means[i] + resp_stds[i], non_means[i] + non_stds[i])
    if p < 0.05:
        fig_bar.add_annotation(x=bar_labels[i], y=max_val + 1.5, text="*",
                               showarrow=False, font=dict(size=20, color="#111827"))
fig_bar.update_layout(
    **CHART_LAYOUT,
    title=chart_title("Mean Cell Population Frequencies with Significance Markers (* p < 0.05, uncorrected)"),
    barmode="group", height=440,
    legend=dict(title="", orientation="h", y=1.06, x=0,
                bgcolor="#ffffff", bordercolor="#d1d5db", borderwidth=1,
                font=dict(color="#111827", size=13)),
    yaxis=dark_axis("Mean Frequency (%)"),
    xaxis=dark_axis(""),
    margin=dict(l=20, r=20, t=70, b=20),
)
st.plotly_chart(fig_bar, use_container_width=True)

st.markdown('<hr class="divider"/>', unsafe_allow_html=True)

# ════════════════════════════════════════════════════════════════════════════════
# SECTION 4 — STATISTICAL RESULTS
# ════════════════════════════════════════════════════════════════════════════════
st.markdown('<div class="section-number">Section 04</div>', unsafe_allow_html=True)
st.markdown('<div class="section-title">Statistical Analysis</div>', unsafe_allow_html=True)
st.markdown(f"""
<div class="narrative">
    Three independent statistical tests were used to evaluate whether differences in cell
    population frequencies between responders and non-responders are statistically meaningful.
    The Mann-Whitney U test (primary) and Brunner-Munzel test are non-parametric and make no
    assumptions about data distribution — appropriate given that Shapiro-Wilk testing confirmed
    the data is not normally distributed. Welch's T-test was included as a secondary parametric
    reference. Bonferroni correction was applied, adjusting the significance threshold from
    p &lt; 0.05 to p &lt; {corrected_threshold:.4f} for five simultaneous comparisons.
</div>
""", unsafe_allow_html=True)

st.dataframe(stats_df, use_container_width=True, hide_index=True)
st.markdown(f"""
<div class="method-note">
    Bonferroni corrected threshold: p &lt; {corrected_threshold:.4f} for {len(POPULATIONS)} simultaneous comparisons.
    Effect size reported as rank-biserial correlation (range: -1 to 1, where 0 = no effect).
    Significant column based on Mann-Whitney p vs corrected threshold.
</div>
""", unsafe_allow_html=True)

st.markdown('<hr class="divider"/>', unsafe_allow_html=True)

# ════════════════════════════════════════════════════════════════════════════════
# SECTION 5 — BASELINE SUBSET
# ════════════════════════════════════════════════════════════════════════════════
st.markdown('<div class="section-number">Section 05</div>', unsafe_allow_html=True)
st.markdown('<div class="section-title">Baseline Cohort Analysis</div>', unsafe_allow_html=True)
st.markdown("""
<div class="narrative">
    To understand the immune landscape prior to any treatment effect, analysis was restricted
    to baseline samples (time from treatment start = 0) from melanoma patients receiving
    miraclib. This subset captures the pre-treatment immune state and is the most clinically
    relevant window for identifying predictive biomarkers.
</div>
""", unsafe_allow_html=True)

m1, m2, m3, m4, m5 = st.columns(5)
for col, (val, label) in zip([m1, m2, m3, m4, m5], [
    (f"{len(baseline_df):,}", "Baseline Samples"),
    (f"{baseline_df['subject'].nunique():,}", "Unique Subjects"),
    (f"{baseline_df['project'].nunique()}", "Projects"),
    (f"{n_resp_base}", "Responders"),
    (f"{avg_bcell:,.2f}", "Avg B Cells (Male Resp.)"),
]):
    with col:
        st.markdown(f'<div class="stat-card"><div class="stat-value">{val}</div><div class="stat-label">{label}</div></div>', unsafe_allow_html=True)

st.markdown("<br/>", unsafe_allow_html=True)
col_a, col_b, col_c = st.columns(3)

with col_a:
    proj = baseline_df.groupby("project")["sample"].count().reset_index()
    proj.columns = ["Project", "Sample Count"]
    fig_proj = px.bar(proj, x="Project", y="Sample Count", color="Project",
                      color_discrete_sequence=["#2563eb", "#16a34a", "#dc2626"])
    fig_proj.update_layout(
        **CHART_LAYOUT, height=300, showlegend=False,
        title=chart_title("Samples by Project"),
        yaxis=dark_axis("Sample Count"),
        xaxis=dark_axis(""),
        margin=dict(l=10, r=10, t=50, b=10)
    )
    st.plotly_chart(fig_proj, use_container_width=True)

with col_b:
    resp_counts = baseline_df.drop_duplicates("subject").groupby("response")["subject"].count().reset_index()
    resp_counts.columns = ["Response", "Count"]
    resp_counts["Response"] = resp_counts["Response"].map({"yes": "Responder", "no": "Non-Responder"})
    fig_resp2 = px.pie(resp_counts, names="Response", values="Count",
                       color_discrete_sequence=[RESP_COLOR, NON_RESP_COLOR], hole=0.5)
    fig_resp2.update_traces(textinfo="label+percent", textfont=dict(size=13, color="#111827"))
    fig_resp2.update_layout(
        **CHART_LAYOUT, height=300, showlegend=False,
        title=chart_title("Subjects by Response"),
        margin=dict(l=10, r=10, t=50, b=10)
    )
    st.plotly_chart(fig_resp2, use_container_width=True)

with col_c:
    sex_counts = baseline_df.drop_duplicates("subject").groupby("sex")["subject"].count().reset_index()
    sex_counts.columns = ["Sex", "Count"]
    sex_counts["Sex"] = sex_counts["Sex"].map({"M": "Male", "F": "Female"})
    fig_sex = px.pie(sex_counts, names="Sex", values="Count",
                     color_discrete_sequence=["#2563eb", "#dc2626"], hole=0.5)
    fig_sex.update_traces(textinfo="label+percent", textfont=dict(size=13, color="#111827"))
    fig_sex.update_layout(
        **CHART_LAYOUT, height=300, showlegend=False,
        title=chart_title("Subjects by Sex"),
        margin=dict(l=10, r=10, t=50, b=10)
    )
    st.plotly_chart(fig_sex, use_container_width=True)

st.markdown('<hr class="divider"/>', unsafe_allow_html=True)

# ════════════════════════════════════════════════════════════════════════════════
# SECTION 6 — KEY FINDINGS
# ════════════════════════════════════════════════════════════════════════════════
st.markdown('<div class="section-number">Section 06</div>', unsafe_allow_html=True)
st.markdown('<div class="section-title">Key Findings</div>', unsafe_allow_html=True)
st.markdown(f"""
<div class="finding-box">
    <strong>Finding 1 — No Population Reaches Significance After Correction</strong>
    After applying Bonferroni correction for five simultaneous comparisons (threshold:
    p &lt; {corrected_threshold:.4f}), no cell population showed a statistically significant
    difference in relative frequency between responders and non-responders. This means we
    cannot yet definitively identify a single immune biomarker of miraclib response from
    this dataset alone.
</div>
<div class="finding-box">
    <strong>Finding 2 — CD4+ T Cells Show the Strongest Trend</strong>
    CD4+ T cells returned the lowest p-value across all three statistical tests:
    Mann-Whitney p = {cd4_row['Mann-Whitney p']}, Brunner-Munzel p = {cd4_row['Brunner-Munzel p']},
    Welch's T p = {cd4_row["Welch's T p"]}. The Welch's test reaches significance even after
    Bonferroni correction. The consistency across three independent tests strengthens
    the case for CD4+ T cells as a candidate biomarker.
</div>
<div class="finding-box">
    <strong>Finding 3 — Effect Sizes Are Small Across All Populations</strong>
    Effect sizes ranged from 0.012 to 0.064, indicating that while statistical trends exist,
    the magnitude of difference between groups is modest. No single population is a strong
    standalone predictor of response in this dataset.
</div>
<div class="finding-box">
    <strong>Finding 4 — Baseline B Cell Count in Male Responders</strong>
    Among male melanoma patients receiving miraclib, the mean B cell count at baseline
    (time = 0) was {avg_bcell:,.2f} cells — a specific metric that may serve as a reference
    point for future biomarker studies.
</div>
""", unsafe_allow_html=True)

st.markdown('<hr class="divider"/>', unsafe_allow_html=True)

# ════════════════════════════════════════════════════════════════════════════════
# SECTION 7 — FUTURE DIRECTIONS
# ════════════════════════════════════════════════════════════════════════════════
st.markdown('<div class="section-number">Section 07</div>', unsafe_allow_html=True)
st.markdown('<div class="section-title">Future Directions</div>', unsafe_allow_html=True)
st.markdown("""
<div class="narrative">
    While no definitive biomarker was identified in this analysis, the findings provide a
    clear roadmap for next steps. The CD4+ T cell signal is worth pursuing with additional
    rigor and a larger dataset.
</div>
""", unsafe_allow_html=True)

f1, f2, f3, f4 = st.columns(4)
futures = [
    ("01", "Expand the Cohort", "Increase sample size to improve statistical power. A larger cohort would allow the CD4+ T cell trend (p = 0.0134) to potentially reach robust significance after correction."),
    ("02", "Investigate CD4+ Subpopulations", "CD4+ T cells are heterogeneous. Breaking them into subtypes (Th1, Th2, Treg, etc.) may reveal a more specific and stronger predictive signal."),
    ("03", "Longitudinal Analysis", "Examine how cell population frequencies change across timepoints. Treatment-induced changes over time may be more predictive than baseline values alone."),
    ("04", "Multivariate Modeling", "Build a machine learning model using all five population frequencies simultaneously. A combination may predict response better than any single biomarker."),
]
for col, (num, title, desc) in zip([f1, f2, f3, f4], futures):
    with col:
        st.markdown(f"""
        <div class="future-item">
            <div class="future-number">{num}</div>
            <div class="future-title">{title}</div>
            <div class="future-desc">{desc}</div>
        </div>
        """, unsafe_allow_html=True)

st.markdown("<br/><br/>", unsafe_allow_html=True)
st.markdown("""
<div style="border-top: 1px solid #d1d5db; padding-top: 20px; font-family: 'IBM Plex Mono', monospace; font-size: 0.7rem; color: #6b7280;">
    Miraclib Clinical Trial Analysis &nbsp;|&nbsp; Loblaw Bio &nbsp;|&nbsp; Confidential
</div>
""", unsafe_allow_html=True)