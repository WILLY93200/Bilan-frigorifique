# app.py
import streamlit as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet
import io, base64

# Constantes
Cp_eau = 4180  # J/kg.K
rho_eau = 1000 # kg/m³

# Conversion °C <-> K
def K(x): return x + 273.15
def C(x): return x - 273.15

# =========================
# FONCTIONS DE CALCUL
# =========================
def calc_cycle(fluide, P_evap_bar, P_cond_bar, 
               T1_C, Tliq_C,
               cond_type, ref_is_outlet,
               T_air_in, T_air_out,
               T_eau_in, T_eau_out,
               eta_isentropic=0.75,
               debit_eau=0.0):

    # Pressions
    P_evap = P_evap_bar * 1e5
    P_cond = P_cond_bar * 1e5

    # Températures
    T1 = K(T1_C)
    Tliq = K(Tliq_C)

    # Saturations
    T_sat_evap = PropsSI("T","P",P_evap,"Q",1,fluide)
    T_sat_cond = PropsSI("T","P",P_cond,"Q",0,fluide)

    # Grandeurs principales
    surchauffe = T1 - T_sat_evap
    sous_refroid = T_sat_cond - Tliq

    # Pincements
    pinc_evap = None
    pinc_cond = None
    if ref_is_outlet:
        pinc_evap = K(T_eau_out) - T_sat_evap
    else:
        pinc_evap = K(T_eau_in) - T_sat_evap

    if cond_type == "Air":
        Tref = K(T_air_out) if ref_is_outlet else K(T_air_in)
        pinc_cond = T_sat_cond - Tref
    else:  # Plaques ou Noyé
        Tref = K(T_eau_out) if ref_is_outlet else K(T_eau_in)
        pinc_cond = T_sat_cond - Tref

    # Cycle
    h1 = PropsSI("H","P",P_evap,"T",T1,fluide)
    s1 = PropsSI("S","P",P_evap,"T",T1,fluide)
    h2s = PropsSI("H","P",P_cond,"S",s1,fluide)
    h2 = h1 + (h2s - h1)/eta_isentropic
    T2 = PropsSI("T","P",P_cond,"H",h2,fluide)
    h3 = PropsSI("H","P",P_cond,"T",Tliq,fluide)
    h4 = h3

    # Débit eau -> puissance utile
    mdot_eau = debit_eau
    Qf = mdot_eau * Cp_eau * (T_eau_in - T_eau_out) / 1000  # kW

    # Débit frigo (approx à partir Qf)
    h_evap_in = PropsSI("H","P",P_evap,"Q",0,fluide)
    h_evap_out = h1
    h_evap_diff = h_evap_out - h_evap_in
    mdot_frigo = Qf*1000 / h_evap_diff if h_evap_diff>0 else 0

    Wc = mdot_frigo*(h2-h1)/1000  # kW
    Qc = Qf + Wc
    COP = Qf/Wc if Wc>0 else None

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

# =========================
# GRAPHIQUES
# =========================
def plot_bars(Qf,Wc,Qc):
    fig,ax=plt.subplots(figsize=(5,3))
    ax.bar(["Qf (utile)","Wc (comp)","Qc (cond)"], [Qf,Wc,Qc], color=["blue","red","green"])
    ax.set_ylabel("Puissance [kW]")
    ax.set_title("Bilan énergétique")
    st.pyplot(fig)

def plot_logph(fluide,res):
    pts=res["points"]
    h=[pts["1"]["h_kJkg"],pts["2"]["h_kJkg"],pts["3"]["h_kJkg"],pts["4"]["h_kJkg"],pts["1"]["h_kJkg"]]
    P=[pts["1"]["P_bar"],pts["2"]["P_bar"],pts["3"]["P_bar"],pts["4"]["P_bar"],pts["1"]["P_bar"]]
    fig,ax=plt.subplots(figsize=(6,4))
    ax.semilogy(h,P,marker="o")
    ax.set_xlabel("h [kJ/kg]")
    ax.set_ylabel("P [bar] (log)")
    ax.set_title(f"log(p)-h cycle – {fluide}")
    ax.grid(True,which="both")
    st.pyplot(fig)

def plot_mollier(fluide,res):
    pts=res["points"]
    # Dôme de saturation
    T_range = np.linspace(PropsSI("Tmin","",0,"",0,fluide)+1, PropsSI("Tcrit","",0,"",0,fluide)-1, 150)
    s_liq,h_liq,s_vap,h_vap=[],[],[],[]
    for T in T_range:
        try:
            s_liq.append(PropsSI("S","T",T,"Q",0,fluide)/1000)
            h_liq.append(PropsSI("H","T",T,"Q",0,fluide)/1000)
            s_vap.append(PropsSI("S","T",T,"Q",1,fluide)/1000)
            h_vap.append(PropsSI("H","T",T,"Q",1,fluide)/1000)
        except: pass

    fig,ax=plt.subplots(figsize=(6,5))
    ax.plot(s_liq,h_liq,'b-',lw=1)
    ax.plot(s_vap,h_vap,'b-',lw=1)
    ax.fill_betweenx(h_liq,s_liq,s_vap,color="lightblue",alpha=0.3)

    # Cycle
    H=[pts["1"]["h_kJkg"],pts["2"]["h_kJkg"],pts["3"]["h_kJkg"],pts["4"]["h_kJkg"],pts["1"]["h_kJkg"]]
    S=[pts["1"]["s_kJkgK"],None,None,None,pts["1"]["s_kJkgK"]]
    try: S[1]=PropsSI("S","P",pts["2"]["P_bar"]*1e5,"H",pts["2"]["h_kJkg"]*1000,fluide)/1000
    except: S[1]=np.nan
    try: S[2]=PropsSI("S","P",pts["3"]["P_bar"]*1e5,"H",pts["3"]["h_kJkg"]*1000,fluide)/1000
    except: S[2]=np.nan
    try: S[3]=PropsSI("S","P",pts["4"]["P_bar"]*1e5,"H",pts["4"]["h_kJkg"]*1000,fluide)/1000
    except: S[3]=np.nan
    ax.plot(S,H,marker="o",color="r",lw=2)

    ax.set_xlabel("s [kJ/kg·K]")
    ax.set_ylabel("h [kJ/kg]")
    ax.set_title(f"Mollier h-s – {fluide}")
    ax.grid(True)
    st.pyplot(fig)

# =========================
# EXPORTS
# =========================
def export_csv(res):
    df=pd.DataFrame([res])
    return df.to_csv(index=False).encode('utf-8')

def export_pdf(res):
    buf=io.BytesIO()
    styles=getSampleStyleSheet()
    doc=SimpleDocTemplate(buf,pagesize=A4)
    story=[Paragraph("Bilan frigorifique",styles["Title"]),Spacer(1,12)]
    for k,v in res.items():
        if isinstance(v,(float,int)):
            txt=f"{k}: {v:.3f}"
        else:
            txt=f"{k}: {v}"
        story.append(Paragraph(txt,styles["Normal"]))
    doc.build(story)
    pdf=buf.getvalue()
    buf.close()
    return pdf

# =========================
# INTERFACE STREAMLIT
# =========================
st.title("📘 Bilan Frigorifique – Application Web")

fluide=st.selectbox("Fluide frigorigène",["R134a","R410A","R32","R407C","R404A","R1234yf","R1234ze","R22","CO2","Ammonia","Propane"])
P_ev=st.number_input("BP (bar)",value=5.94)
P_co=st.number_input("HP (bar)",value=24.62)
T1=st.number_input("T sortie évap (°C)",value=10.0)
Tliq=st.number_input("T ligne liquide (°C)",value=35.0)

cond_type=st.selectbox("Condenseur",["Air","Plaques (eau)","Noyé"])
ref_is_outlet=st.radio("Référence pincement",["Entrée","Sortie"])=="Sortie"

Tair_in=st.number_input("T air entrée (°C)",value=25.0)
Tair_out=st.number_input("T air sortie (°C)",value=32.0)
Teau_in=st.number_input("T eau entrée (°C)",value=12.0)
Teau_out=st.number_input("T eau sortie (°C)",value=7.0)

eta=st.slider("Rendement isentropique compresseur",0.5,0.9,0.75,0.01)

unit=st.radio("Unité de débit eau",["m³/h","L/s"])
flow=st.number_input("Débit d'eau",value=2.0)

if st.button("Calculer"):
    if unit=="m³/h":
        mdot=flow/3.6*rho_eau
    else:
        mdot=flow*rho_eau

    res=calc_cycle(fluide,P_ev,P_co,T1,Tliq,cond_type,ref_is_outlet,Tair_in,Tair_out,Teau_in,Teau_out,eta,mdot)

    st.subheader("Résultats")
    st.write(f"**Surchauffe**: {res['surchauffe_K']:.2f} K")
    st.write(f"**Sous-refroidissement**: {res['sous_refroid_K']:.2f} K")
    st.write(f"**Pincement évap**: {res['pincement_evap_K']:.2f} K")
    st.write(f"**Pincement condenseur**: {res['pincement_cond_K']:.2f} K")
    st.write(f"**Qf**: {res['Qf_kW']:.2f} kW, **Wc**: {res['Wc_kW']:.2f} kW, **Qc**: {res['Qc_kW']:.2f} kW, **COP**: {res['COP']:.2f}")

    plot_bars(res["Qf_kW"],res["Wc_kW"],res["Qc_kW"])
    plot_logph(fluide,res)
    plot_mollier(fluide,res)

    # Exports
    csv=export_csv(res)
    pdf=export_pdf(res)
    st.download_button("⬇️ Télécharger CSV",csv,"resultats.csv","text/csv")
    st.download_button("⬇️ Télécharger PDF",pdf,"resultats.pdf","application/pdf")
