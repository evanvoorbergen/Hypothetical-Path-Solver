
import pandas as pd
from scipy.integrate import quad
import re
import numpy as np


Table_B1 = pd.read_excel('..\data\IPT_Tabel_B1.xlsx').fillna(0)
Table_B1 = Table_B1.apply(lambda col: col.str.strip() if col.dtype == "object" else col) #remove spaces

# make numerical columns float
num_cols_B1 = ["Molaire massa [g/mol]", "SG (20°/4°)", "Tm [°C]", "Hm(Tm) [kJ/mol]", "Tb [°C]", "Hv(Tb) [kJ/mol]", "Tc [K]", "Pc [atm]"]
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


Table_B1 

Table_B2 = pd.read_excel('..\data\IPT_Tabel_B2.xlsx')
Table_B2 = Table_B2.apply(lambda col: col.str.strip()) #remove spaces
Table_B2 = Table_B2.fillna(0) # fill NaN w 0

num_cols_B2 = ['Molaire massa', 'Vorm', 'a (•10^3)', 'b (•10^5)', 'c (•10^8)', 'd (•10^12)']

Table_B2[num_cols_B2] = Table_B2[num_cols_B2].apply(pd.to_numeric, errors="coerce")

Table_B2

class Molecule:
    def __init__(self, name, T, P, F, Table_B1, Table_B2):
        self.name = name
        
        self.T = T
        self.P = P
        self.F = F

        self.Table_B1 = Table_B1[(Table_B1["Stofnaam (NL)"] == self.name) | (Table_B1["Stofnaam (EN)"] == self.name) | (Table_B1["Formule"] == self.name)].iloc[0]       # Selects row of Table B1 for the correct compound
        self.Table_B2 = Table_B2[((Table_B2["Stofnaam (NL)"] == self.name) | (Table_B2["Stofnaam (EN)"] == self.name) | (Table_B2["Formule"] == self.name)) & (Table_B2["Staat"] == self.F)].iloc[0]         # Selects row of Table B2 for the correct compound and phase
    

    ### Temperature Change

    # Get Cp (not integrated yet!)
    def Cp(self, T = None):
        if T == None:
            T = self.T

        coeff = self.Table_B2
        a = coeff['a (•10^3)']*1e-3
        b = coeff['b (•10^5)']*1e-5
        c = coeff['c (•10^8)']*1e-8
        d = coeff['d (•10^12)']*1e-12
        
        return  a + b*T + c*T**2 + d*T**3
    
    # Integrate CpdT between T boundaries
    def temp_change(self, new_T):
        dH = quad(lambda T: self.Cp(T), self.T, new_T)[0] #kJ/mol

        # Adjust to New Temperature
        self.T = new_T
        return dH
    

    ### Phase Change

    def phase_change(self, new_F):
        solid = 's' or 'c'      # idk why but solid is given as 'c'

        # Melting and Freezing
        if self.F == solid and new_F == 'l':
            dH = self.Table_B1['Hm(Tm) [kJ/mol]']
        elif self.F == 'l' and new_F == solid:
            dH = -self.Table_B1['Hm(Tm) [kJ/mol]']
        # Vaporization and Condensation
        elif self.F == 'l' and new_F == 'g':
            dH = self.Table_B1['Hv(Tb) [kJ/mol]']
        elif self.F == 'g' and new_F == 'l':
            dH = -self.Table_B1['Hv(Tb) [kJ/mol]']
        # No Phase Change
        elif self.F == new_F:
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

        
# Try it Out
name = 'H2O'
F_start = 'g'
F_end = 'g'
T_start = 25 # C
T_end = 546 # C
P_start = 1  # atm
P_end = 1    # atm

# Create molecule
water = Molecule(name, T_start, P_start, F_start, Table_B1, Table_B2)
path = Path(water, T_end, P_end, F_end)
path.build_path()
print(path.total_enthalpy()) 


### STREAMLIT ###

import streamlit as st

st.title("Thermodynamic Path Builder")

all_molecules = np.array(Table_B2["Formule"])

name = st.selectbox("Molecule", all_molecules)
T_start = st.number_input("Start Temperature (°C)", value=25)
T_end = st.number_input("End Temperature (°C)", value=100)

F_start = st.selectbox("Start Phase", ["s", "l", "g"])
F_end = st.selectbox("End Phase", ["s", "l", "g"])

if st.button("Calculate Path"):
    mol = Molecule(name, T_start, 1, F_start, Table_B1, Table_B2)
    path = Path(mol, T_end, 1, F_end)
    path.build_path()

    st.write("Total ΔH:", path.total_enthalpy())

    for step in path.steps:
        st.write(step.dH)