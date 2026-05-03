import re
import pandas as pd
from scipy.integrate import quad
import numpy as np


Table_B1 = pd.read_excel(r'../data/IPT_Tabel_B1.xlsx').fillna(0)
Table_B1 = Table_B1.apply(
    lambda col: col.str.strip() if col.dtype == "object" else col
    ) #remove spaces

# make numerical columns float
num_cols_B1 = ["Molaire massa [g/mol]", 
               "SG (20°/4°)", 
               "Tm [°C]", 
               "Hm(Tm) [kJ/mol]", 
               "Tb [°C]", 
               "Hv(Tb) [kJ/mol]", 
               "Tc [K]", 
               "Pc [atm]"
               ]
Table_B1[num_cols_B1] = Table_B1[num_cols_B1].apply(pd.to_numeric, errors="coerce")

# make Hc string
Table_B1["Hc° [kJ/mol]"] = Table_B1["Hc° [kJ/mol]"].replace("-", np.nan)

# Extra Pre-Processing for Hf and Hc

# Split Heat of Formation Values for Different Phases
for row in range(len(Table_B1)):
    string = Table_B1["Hf° [kJ/mol]"].iloc[row]
    string = string.replace('`', '')
    splits = string.split()
    for i in range(len(splits)):
        s = splits[i]
        s = re.match(r"(-?\d+\.?\d*)\((\w+)\)", s)
        if s:
            value = float(s.group(1))
            phase = str(s.group(2))
            Table_B1.loc[row, f"Hf° ({phase}) [kJ/mol]"] = value


# Split Heat of Combustion Values for Different Phases
for row in range(len(Table_B1)):
    string = Table_B1["Hc° [kJ/mol]"].iloc[row]

    if pd.isna(string):
        continue

    string = string.replace('`', '')
    splits = string.split()
    for i in range(len(splits)):
        s = splits[i]
        s = re.match(r"(-?\d+\.?\d*)\((\w+)\)", s)
        if s:
            value = float(s.group(1))
            phase = str(s.group(2))
            Table_B1.loc[row, f"Hc° ({phase}) [kJ/mol]"] = value

# Drop original Hc and Hf columns
Table_B1 = Table_B1.drop(['Hf° [kJ/mol]', 'Hc° [kJ/mol]'], axis=1)


Table_B2 = pd.read_excel(r"../data/IPT_Tabel_B2.xlsx")
Table_B2 = Table_B2.apply(lambda col: col.str.strip()) #remove spaces # if col.dtype == str else col
Table_B2 = Table_B2.fillna(0) # fill NaN w 0

num_cols_B2 = ['Molaire massa', 'Vorm', 'a (•10^3)', 'b (•10^5)', 'c (•10^8)', 'd (•10^12)']

Table_B2[num_cols_B2] = Table_B2[num_cols_B2].apply(pd.to_numeric, errors="coerce")

Table_B2 = Table_B2.fillna(0) # fill NaN w 0


class Molecule:
    def __init__(self, name, T, P, F, Table_B1, Table_B2):
        self.name = name
        
        self.T = T
        self.P = P

        if F == 's':
            self.F = 'c'
        else:
            self.F = F

        self.Table_B1 = Table_B1[(Table_B1["Stofnaam (NL)"] == self.name) | (Table_B1["Stofnaam (EN)"] == self.name) | (Table_B1["Formule"] == self.name)].iloc[0]       # Selects row of Table B1 for the correct compound
        try:
            self.Table_B2 = Table_B2[((Table_B2["Stofnaam (NL)"] == self.name) | (Table_B2["Stofnaam (EN)"] == self.name) | (Table_B2["Formule"] == self.name)) & (Table_B2["Staat"] == self.F)].iloc[0]         # Selects row of Table B2 for the correct compound and phase
        except:
            raise ValueError("There is no defined heat capacity function for this molecule.")

    ### Temperature Change
    
    def temp_change(self, new_T):

        coeff = self.Table_B2

        a = coeff['a (•10^3)'] * 1e-3
        b = coeff['b (•10^5)'] * 1e-5
        c = coeff['c (•10^8)'] * 1e-8
        d = coeff['d (•10^12)'] * 1e-12

        T1 = self.T
        T2 = new_T

        # Convert Temp to K if Cp is defined in K
        if coeff['Eenheid van temperatuur'] == 'K':
            T1 = self.T + 273.15 #K
            T2 = new_T + 273.15 #K


        if coeff['Vorm'] == 1:
            def H(T):
                return a*T + (b/2)*T**2 + (c/3)*T**3 + (d/4)*T**4   # Integrated form of a + b*T + c*T**2 + d*T**3
        elif coeff['Vorm'] == 2:
            def H(T):
                return a*T + (b/2)*T**2 - c/T                       # Integrated form of a + b*T + c*T**(-2)
        else:
            raise ValueError("Cp function is not of Form 1 or Form 2")

        # Calculate enthalpy change
        dH = (H(T2) - H(T1))  # kJ/mol

        # Update T (Celsius!)
        self.T = new_T
        return dH


    ### Phase Change

    def phase_change(self, new_F):
        solid = ['s', 'c']      # idk why but solid is given as 'c'

        # Melting and Freezing
        if self.F in solid and new_F == 'l':
            dH = self.Table_B1['Hm(Tm) [kJ/mol]']
        elif self.F == 'l' and new_F in solid:
            dH = -self.Table_B1['Hm(Tm) [kJ/mol]']
        # Vaporization and Condensation
        elif self.F == 'l' and new_F == 'g':
            dH = self.Table_B1['Hv(Tb) [kJ/mol]']
        elif self.F == 'g' and new_F == 'l':
            dH = -self.Table_B1['Hv(Tb) [kJ/mol]']
        # No Phase Change
        elif self.F == new_F or (self.F in solid and new_F in solid):
            dH = 0
        # Error for Other Changes
        else:
            raise ValueError('There is no enthalpy provided for this phase transition.')

        # Adjust to new phase
        self.F = new_F
        return dH
    

class Step:
    def __init__(self, molecule, new_T = None, new_F = None, new_P = None):
        self.molecule = molecule

        self.current_T = molecule.T
        self.current_P = molecule.P
        self.current_F = molecule.F

        self.new_T = new_T or molecule.T
        self.new_P = new_P or molecule.P
        self.new_F = new_F or molecule.F

        self.dH = 0

        if self.current_T != self.new_T:
            self.dH += molecule.temp_change(new_T)
        
        if self.current_F != self.new_F:
            self.dH += molecule.phase_change(new_F)
        
class Path:
    def __init__(self, molecule, end_T, end_P, end_F):
        self.steps = []

        self.molecule = molecule

        self.end_T = end_T
        self.end_P = end_P
        self.end_F = end_F

    def add_step(self, step):
        self.steps.append(step)

    def total_enthalpy(self):
        enthalpy = 0
        for step in self.steps:
            enthalpy += step.dH
        return enthalpy

    def build_path(self):
        # Go to Tmelt and Phase Change to Liquid
        if self.molecule.F in ['c', 's'] and self.end_F in ['l', 'g']:
            self.add_step(Step(self.molecule, new_T = self.molecule.Table_B1['Tm [°C]']))
            self.add_step(Step(self.molecule, new_F = 'l'))

        # Go to Tmelt and Phase Change to Solid
        if self.molecule.F == 'l' and self.end_F in ['c', 's']:
            self.add_step(Step(self.molecule, new_T = self.molecule.Table_B1['Tm [°C]']))
            self.add_step(Step(self.molecule, new_F = 'c'))

        # Go to Tvap and Phase Change to Liquid if Ending in Solid
        if self.molecule.F == 'g' and self.end_F in ['c', 's']:
            self.add_step(Step(self.molecule, new_T = self.molecule.Table_B1['Tb [°C]']))
            self.add_step(Step(self.molecule, new_F = 'l'))

        # Go to Tvap and Change Phase
        if {self.molecule.F, self.end_F} == {'l', 'g'}:
            self.add_step(Step(self.molecule, new_T = self.molecule.Table_B1['Tb [°C]']))
            if self.molecule.F == 'l' and self.end_F == 'g':
                self.add_step(Step(self.molecule, new_F='g'))
            elif self.molecule.F == 'g' and self.end_F == 'l':
                self.add_step(Step(self.molecule, new_F='l'))
    
        # Final Temp Change
        if self.molecule.F == self.end_F and self.molecule.T != self.end_T:
            self.add_step(Step(self.molecule, new_T = self.end_T))



### STREAMLIT ###

# import streamlit as st

# st.title("Thermodynamic Path Builder")

# all_molecules = np.array(Table_B2["Formule"])

# name = st.selectbox("Molecule", all_molecules)
# T_start = st.number_input("Start Temperature (°C)", value=25)
# T_end = st.number_input("End Temperature (°C)", value=100)

# F_start = st.selectbox("Start Phase", ["s", "l", "g"])
# F_end = st.selectbox("End Phase", ["s", "l", "g"])

# if st.button("Calculate Path"):
#     mol = Molecule(name, T_start, 1, F_start, Table_B1, Table_B2)
#     path = Path(mol, T_end, 1, F_end)
#     path.build_path()

#     st.write("Total ΔH:", path.total_enthalpy())

#     for step in path.steps:
#         st.write(step.dH)
        
### STREAMLIT ###

import streamlit as st
import streamlit.components.v1 as components

st.title("Thermodynamic Path Builder")

all_molecules = np.array(Table_B2["Formule"])

name = st.selectbox("Molecule", all_molecules)
T_start = st.number_input("Start Temperature (°C)", value=25)
T_end = st.number_input("End Temperature (°C)", value=100)

F_start = st.selectbox("Start Phase", ["s", "l", "g"])
F_end = st.selectbox("End Phase", ["s", "l", "g"])

PHASE_LABELS = {"s": "solid", "c": "solid", "l": "liquid", "g": "gas"}

def get_step_label(step):
    """Determine the thermodynamic label for a step's arrow."""
    T_changed = step.current_T != step.new_T
    F_changed = step.current_F != step.new_F

    if T_changed and not F_changed:
        phase = step.current_F
        if phase in ("s", "c"):
            return "Cp,s·dT"
        elif phase == "l":
            return "Cp,l·dT"
        elif phase == "g":
            return "Cp,g·dT"

    if F_changed and not T_changed:
        solid = ("s", "c")
        if step.current_F in solid and step.new_F == "l":
            return "ΔH_fus"
        elif step.current_F == "l" and step.new_F in solid:
            return "−ΔH_fus"
        elif step.current_F == "l" and step.new_F == "g":
            return "ΔH_vap"
        elif step.current_F == "g" and step.new_F == "l":
            return "−ΔH_vap"

    return "ΔH"

def build_diagram_html(molecule_name, steps):
    BOX_W, BOX_H = 160, 80
    ARROW_W = 120
    STEP_W = BOX_W + ARROW_W
    TOTAL_BOXES = len(steps) + 1
    SVG_W = TOTAL_BOXES * BOX_W + len(steps) * ARROW_W + 40
    SVG_H = 200

    # Collect all states (box before each step + final box)
    states = []
    for step in steps:
        states.append((step.current_T, step.current_F))
    last = steps[-1]
    states.append((last.new_T, last.new_F))

    def phase_color(f):
        return {
            "s": ("#dbeafe", "#1e40af"),   # blue: solid
            "c": ("#dbeafe", "#1e40af"),
            "l": ("#d1fae5", "#065f46"),   # green: liquid
            "g": ("#fef3c7", "#92400e"),   # amber: gas
        }.get(f, ("#f3f4f6", "#374151"))

    svg_parts = [
        f'<svg width="{SVG_W}" height="{SVG_H}" viewBox="0 0 {SVG_W} {SVG_H}" '
        f'xmlns="http://www.w3.org/2000/svg" font-family="sans-serif">',
        '<defs><marker id="arr" viewBox="0 0 10 10" refX="8" refY="5" '
        'markerWidth="6" markerHeight="6" orient="auto">'
        '<path d="M2 1L8 5L2 9" fill="none" stroke="#6b7280" '
        'stroke-width="1.5" stroke-linecap="round" stroke-linejoin="round"/>'
        '</marker></defs>',
    ]

    CY = SVG_H // 2  # vertical center
    x = 20

    for i, (T, F) in enumerate(states):
        fill, text_col = phase_color(F)
        phase_str = PHASE_LABELS.get(F, F)

        # Box
        svg_parts.append(
            f'<rect x="{x}" y="{CY - BOX_H//2}" width="{BOX_W}" height="{BOX_H}" '
            f'rx="10" fill="{fill}" stroke="{text_col}" stroke-width="1.2"/>'
        )
        # Molecule name
        svg_parts.append(
            f'<text x="{x + BOX_W//2}" y="{CY - 14}" text-anchor="middle" '
            f'font-size="15" font-weight="600" fill="{text_col}">{molecule_name}</text>'
        )
        # Temperature
        svg_parts.append(
            f'<text x="{x + BOX_W//2}" y="{CY + 6}" text-anchor="middle" '
            f'font-size="13" fill="{text_col}">T = {T:.1f} °C</text>'
        )
        # Phase
        svg_parts.append(
            f'<text x="{x + BOX_W//2}" y="{CY + 22}" text-anchor="middle" '
            f'font-size="12" fill="{text_col}">({phase_str})</text>'
        )

        # Arrow to next box (if not last state)
        if i < len(steps):
            step = steps[i]
            label = get_step_label(step)
            dH_val = step.dH
            sign = "+" if dH_val >= 0 else ""
            ax1 = x + BOX_W
            ax2 = ax1 + ARROW_W

            # Arrow line
            svg_parts.append(
                f'<line x1="{ax1}" y1="{CY}" x2="{ax2 - 6}" y2="{CY}" '
                f'stroke="#6b7280" stroke-width="1.5" marker-end="url(#arr)"/>'
            )
            # Step type label (above arrow)
            svg_parts.append(
                f'<text x="{(ax1 + ax2)//2}" y="{CY - 14}" text-anchor="middle" '
                f'font-size="12" font-weight="600" fill="#374151">{label}</text>'
            )
            # ΔH value (below arrow)
            svg_parts.append(
                f'<text x="{(ax1 + ax2)//2}" y="{CY + 18}" text-anchor="middle" '
                f'font-size="11" fill="#6b7280">{sign}{dH_val:.2f} kJ/mol</text>'
            )

            x = ax2

    svg_parts.append('</svg>')
    svg_html = "\n".join(svg_parts)

    # Wrap in scrollable div in case diagram is wide
    return f"""
    <div style="overflow-x: auto; padding: 12px 0;">
      {svg_html}
    </div>
    """

if st.button("Calculate Path"):
    try:
        mol = Molecule(name, T_start, 1, F_start, Table_B1, Table_B2)
        path = Path(mol, T_end, 1, F_end)
        path.build_path()

        st.write(f"**Total ΔH:** {path.total_enthalpy():.4f} kJ/mol")

        if path.steps:
            html = build_diagram_html(name, path.steps)
            components.html(html, height=220, scrolling=True)

            # Also show breakdown table
            st.subheader("Step breakdown")
            for i, step in enumerate(path.steps):
                label = get_step_label(step)
                phase_before = PHASE_LABELS.get(step.current_F, step.current_F)
                phase_after  = PHASE_LABELS.get(step.new_F, step.new_F)
                st.write(
                    f"**Step {i+1}** ({label}): "
                    f"{step.current_T:.1f}°C {phase_before} → "
                    f"{step.new_T:.1f}°C {phase_after} | "
                    f"ΔH = {step.dH:.4f} kJ/mol"
                )
        else:
            st.info("No steps generated — start and end states may be identical.")

    except ValueError as e:
        st.error(str(e))
