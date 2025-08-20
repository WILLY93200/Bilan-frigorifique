# app.py – Bilan frigorifique (version complète avec fluides récents)

import streamlit as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet
import io

# ---------------------------
# Constantes & utilitaires
# ---------------------------
Cp_eau = 4180   # J/kg.K
rho_eau = 1000  # kg/m3

def K(x): return x + 273.15
def C(x): return x - 273.15

def fmt(x, digits=2, unit=""):
    """Formatte proprement une valeur éventuellement None/NaN/Inf."""
    try:
        if x is None or (isinstance(x, float) and (np.isnan(x) or np.isinf(x))):
            return "—"
        return f"{x:.{digits}f}{unit}"
    except Exception:
        return "—"

# --- Liste élargie de fluides CoolProp (libellé -> clé CoolProp) ---
FLUIDS = {
    # HFC classiques
    "R22": "R22",
    "R134a": "R134a",
    "R404A": "R404A",
    "R407A": "R407A",
    "R407C": "R407C",
    "R407F": "R407F",
    "R410A": "R410A",
    "R507A": "R507A",

    # HFO et mélanges récents bas GWP
    "R1234yf": "R1234yf",
    "R1234ze(E)": "R1234ze",
    "R1233zd(E)": "R1233zd(E)",
    "R450A (XP10)": "R450A",
    "R513A (XP10 A)": "R513A",
    "R448A (Solstice N40)": "R448A",
    "R449A (Opteon XP40)": "R449A",
    "R452A": "R452A",
    "R452B": "R452B",
    "R454A": "R454A",
    "R454B": "R454B",
    "R454C": "R454C",
    "R455A": "R455A",

    # Naturels
    "CO2 (R744)": "CO2",
    "Ammonia (R717)": "Ammonia",
    "Propane (R290)": "Propane",
    "Isobutane (R600a)": "IsoButane",
    "n-Butane (R600)": "n-Butane",
    "Propylene (R1270)": "Propylene",
}

def fluid_supported(coolprop_key: str) -> bool:
    """Vérifie rapidement si CoolProp connaît ce fluide dans l'environnement courant."""
    try:
        PropsSI("Tcrit", "", 0, "", 0, coolprop_key)
        return True
    except Exception:
        return False


# ---------------------------
# Calculs thermodynamiques
# ---------------------------
def calc_cycle(fluide, P_evap_bar, P_cond_bar,
               T1_C, Tliq_C,
               cond_type, ref_is_outlet,
               T_air_in, T_air_out,
               T_eau_in, T_eau_out,
               eta_isentropic=0.75,
               debit_eau=0.0):

    # Pressions (Pa)
    P_evap = P_evap_bar * 1e5
    P_cond = P_cond_bar * 1e5

    # Températures (K)
    T1 = K(T1_C)       # point 1 : sortie évap (gaz)
    Tliq = K(Tliq_C)   # point 3 : liquide sous-refroidi

    # Saturations
    T_sat_evap = PropsSI("T", "P", P_evap, "Q", 1, fluide)
    T_sat_cond = PropsSI("T", "P", P_cond, "Q", 0, fluide)

    # Grandeurs de base
    surchauffe   = T1 - T_sat_evap
    sous_refroid = T_sat_cond - Tliq

    # Pincements
    pinc_evap = (K(T_eau_out) if ref_is_outlet else K(T_eau_in)) - T_sat_evap

    if cond_type == "Air":
        Tref = K(T_air_out) if ref_is_outlet else K(T_air_in)
        pinc_cond = T_sat_cond - Tref
    else:  # Plaques (eau) ou Noyé
        Tref = K(T_eau_out) if ref_is_outlet else K(T_eau_in)
        pinc_cond = T_sat_cond - Tref

    # Points cycle
    h1 = PropsSI("H", "P", P_evap, "T", T1, fluide)
    s1 = PropsSI("S", "P", P_evap, "T", T1, fluide)
    h2s = PropsSI("H", "P", P_cond, "S", s1, fluide)
    h2  = h1 + (h2s - h1) / max(eta_isentropic, 1e-6)
    T2  = PropsSI("T", "P", P_cond, "H", h2, fluide)
    h3  = PropsSI("H", "P", P_cond, "T", Tliq, fluide)
    h4  = h3  # détente iso-enthalpique

    # Puissance utile via débit d'eau (kg/s)
    mdot_eau = debit_eau
    Qf = mdot_eau * Cp_eau * (T_eau_in - T_eau_out) / 1000  # kW

    # Estimation du débit frigo via l'évaporateur
    try:
        h_evap_in  = PropsSI("H", "P", P_evap, "Q", 0, fluide)
        h_evap_out = h1
        h_evap_diff = max(h_evap_out - h_evap_in, 0.0)
        mdot_frigo  = (Qf * 1000 / h_evap_diff) if h_evap_diff > 0 else 0.0
    except Exception:
        mdot_frigo = 0.0

    Wc = mdot_frigo * max(h2 - h1, 0.0) / 1000  # kW
    Qc = Qf + Wc
    COP = (Qf / Wc) if Wc > 0 else None

    return {
        "fluide": fluide,
        "T_sat_evap_C": C(T_sat_evap),
        "T_sat_cond_C": C(T_sat_cond),
        "surchauffe_K": surchauffe,
        "sous_refroid_K": sous_refroid,
        "pincement_evap_K": pinc_evap,
        "pincement_cond_K": pinc_cond,
        "Qf_kW": Qf, "Wc_kW": Wc, "Qc_kW": Qc, "COP": COP,
        "points": {
            "1": {"P_bar": P_evap_bar, "T_C": C(T1), "h_kJkg": h1/1000, "s_kJkgK": s1/1000},
            "2": {"P_bar": P_cond_bar, "T_C": C(T2), "h_kJkg": h2/1000},
            "3": {"P_bar": P_cond_bar, "T_C": C(Tliq), "h_kJkg": h3/1000},
            "4": {"P_bar": P_evap_bar, "h_kJkg": h4/1000}
        }
    }


# ---------------------------
# Graphiques
# ---------------------------
def plot_bars(Qf, Wc, Qc):
    fig, ax = plt.subplots(figsize=(5, 3))
    ax.bar(["Qf (utile)", "Wc (comp)", "Qc (cond)"], [Qf, Wc, Qc], color=["blue", "red", "green"])
    ax.set_ylabel("Puissance [kW]")
    ax.set_title("Bilan énergétique")
    st.pyplot(fig)

def plot_logph(fluide, res):
    pts = res["points"]
    h = [pts["1"]["h_kJkg"], pts["2"]["h_kJkg"], pts["3"]["h_kJkg"], pts["4"]["h_kJkg"], pts["1"]["h_kJkg"]]
    P = [pts["1"]["P_bar"],   pts["2"]["P_bar"],   pts["3"]["P_bar"],   pts["4"]["P_bar"],   pts["1"]["P_bar"]]
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.semilogy(h, P, marker="o")
    ax.set_xlabel("h [kJ/kg]")
    ax.set_ylabel("P [bar] (log)")
    ax.set_title(f"log(p)-h cycle – {fluide}")
    ax.grid(True, which="both")
    st.pyplot(fig)

def plot_mollier(fluide, res):
    pts = res["points"]
    T_range = np.linspace(PropsSI("Tmin","",0,"",0,fluide)+1,
                          PropsSI("Tcrit","",0,"",0,fluide)-1, 150)
    s_liq, h_liq, s_vap, h_vap = [], [], [], []
    for T in T_range:
        try:
            s_liq.append(PropsSI("S","T",T,"Q",0,fluide)/1000)
            h_liq.append(PropsSI("H","T",T,"Q",0,fluide)/1000)
            s_vap.append(PropsSI("S","T",T,"Q",1,fluide)/1000)
            h_vap.append(PropsSI("H","T",T,"Q",1,fluide)/1000)
        except Exception:
            pass

    fig, ax = plt.subplots(figsize=(6, 5))
    if s_liq and s_vap:
        ax.plot(s_liq, h_liq, 'b-', lw=1)
        ax.plot(s_vap, h_vap, 'b-', lw=1)
        ax.fill_betweenx(h_liq, s_liq, s_vap, color="lightblue", alpha=0.3)

    H = [pts["1"]["h_kJkg"], pts["2"]["h_kJkg"], pts["3"]["h_kJkg"], pts["4"]["h_kJkg"], pts["1"]["h_kJkg"]]
    S = [pts["1"]["s_kJkgK"], None, None, None, pts["1"]["s_kJkgK"]]
    try: S[1] = PropsSI("S","P",pts["2"]["P_bar"]*1e5,"H",pts["2"]["h_kJkg"]*1000,fluide)/1000
    except Exception: S[1] = np.nan
    try: S[2] = PropsSI("S","P",pts["3"]["P_bar"]*1e5,"H",pts["3"]["h_kJkg"]*1000,fluide)/1000
    except Exception: S[2] = np.nan
    try: S[3] = PropsSI("S","P",pts["4"]["P_bar"]*1e5,"H",pts["4"]["h_kJkg"]*1000,fluide)/1000
    except Exception: S[3] = np.nan

    ax.plot(S, H, marker="o", color="r", lw=2)
    ax.set_xlabel("s [kJ/kg·K]")
    ax.set_ylabel("h [kJ/kg]")
    ax.set_title(f"Mollier h-s – {fluide}")
    ax.grid(True)
    st.pyplot(fig)


# ---------------------------
# Exports
# ---------------------------
def export_csv(res):
    df = pd.DataFrame([res])
    return df.to_csv(index=False).encode("utf-8")

def export_pdf(res):
    buf = io.BytesIO()
    styles = getSampleStyleSheet()
    doc = SimpleDocTemplate(buf, pagesize=A4)
    story = [Paragraph("Bilan frigorifique", styles["Title"]), Spacer(1, 12)]
    for k, v in res.items():
        if isinstance(v, (float, int)):
            txt = f"{k}: {v:.3f}"
        else:
            txt = f"{k}: {v}"
        story.append(Paragraph(txt, styles["Normal"]))
    doc.build(story)
    pdf = buf.getvalue()
    buf.close()
    return pdf


# ---------------------------
# Interface Streamlit
# ---------------------------
st.title("📘 Bilan Frigorifique – Application Web")

labels = list(FLUIDS.keys())
default_idx = labels.index("R410A") if "R410A" in labels else 0
label = st.selectbox("Fluide frigorigène", labels, index=default_idx)
fluide = FLUIDS[label]
if not fluid_supported(fluide):
    st.error(f"Le fluide « {label} » ({fluide}) n’est pas supporté par CoolProp ici.")
    st.stop()

P_ev = st.number_input("BP (bar)", value=5.94)
P_co = st.number_input("HP (bar)", value=24.62)
T1   = st.number_input("T sortie évap (°C)", value=10.0)
Tliq = st.number_input("T ligne liquide (°C)", value=35.0)

cond_type     = st.selectbox("Condenseur", ["Air", "Plaques (eau)", "Noyé"])
ref_is_outlet = (st.radio("Référence pincement", ["Entrée", "Sortie"]) == "Sortie")

Tair_in  = st.number_input("T air entrée (°C)", value=25.0)
Tair_out = st.number_input("T air sortie (°C)", value=32.0)
Teau_in  = st.number_input("T eau entrée (°C)", value=12.0)
Teau_out = st.number_input("T eau sortie (°C)", value=7.0)

eta  = st.slider("Rendement isentropique compresseur", 0.5, 0.9, 0.75, 0.01)
unit = st.radio("Unité de débit eau", ["m³/h", "L/s"])
flow = st.number_input("Débit d'eau", value=2.0)

if st.button("Calculer"):
    if unit == "m³/h":
        mdot = flow/3.6 * rho_eau
    else:
        mdot = flow * rho_eau

    if Teau_in <= Teau_out:
        st.warning("⚠️ Pour un évaporateur, **T eau entrée** doit être > **T eau sortie** (ΔT > 0).")

    res = calc_cycle(fluide, P_ev, P_co, T1, Tliq,
                     cond_type, ref_is_outlet,
                     Tair_in, Tair_out, Teau_in, Teau_out,
                     eta, mdot)

    st.subheader("Résultats")
    st.write(f"**Surchauffe** : {fmt(res['surchauffe_K'])} K")
    st.write(f"**Sous-refroidissement** : {fmt(res['sous_refroid_K'])} K")
    st.write(f"**Pincement évap** : {fmt(res['pincement_evap_K'])} K")
    st.write(f"**Pincement condenseur** : {fmt(res['pincement_cond_K'])} K")
    st.write(
        f"**Qf** : {fmt(res['Qf_kW'])} kW  |  "
        f"**Wc** : {fmt(res['Wc_kW'])} kW  |  "
        f"**Qc** : {fmt(res['Qc_kW'])} kW  |  "
        f"**COP** : {fmt(res['COP'])}"
    )

    plot_bars(res.get("Qf_kW", 0) or 0, res.get("Wc_kW", 0) or 0, res.get("Qc_kW", 0) or 0)
    plot_logph(fluide, res)
    plot_mollier(fluide, res)

    csv = export_csv(res)
    pdf = export_pdf(res)
    st.download_button("⬇️ Télécharger CSV", csv, "resultats.csv", "text/csv")
    st.download_button("⬇️ Télécharger PDF", pdf, "resultats.pdf", "application/pdf")
