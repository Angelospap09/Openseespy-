# Openseespy Timoshenko 


import openseespy.opensees as ops
import math
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import opsvis as opsv
from vfo import vfo
from streng.codes.eurocodes.ec8.cls.seismic_action.spectra import SpectraEc8
from streng.ppp.sections.geometry.tee import TeeSectionGeometry
from streng.ppp.sections.geometry.rectangular import RectangularSectionGeometry
from streng.codes.eurocodes.ec2.raw.ch5.geometric_data import effective_width


# Δεδομένα      
# Μονάδες: kN, m, C

# Διαστάσεις

L = 7             # Μήκος δοκού (m)
H = 3             # Mήκος υποστυλώματος (m)
a = L/4           # Απόσταση φορτίου από ακραίο υποστύλωμα (m)
L1 = L/2          # Μήκος τμήματος δοκού μεταξύ υποστυλωμάτων (m)
cnom = 0.05       # Επικάλυψη (m)

# Δοκός - Τυπική διατομή Τ
bw = 0.3         
hf = 0.15        
hb = 0.5         

# Για μονό άνοιγμα με πλαισιακή σύνδεση (μονολιθική με υποστυλώματα)
l0 = effective_width.l0(l2=L, zero_moments_case=1)  
b1 = b2 = L/2 - bw/2  
beff1 = beff2 = effective_width.beffi(b1, l0)
_b = bw + b1 + b2
beff = effective_width.beff(bw, beff1, beff2, _b)

# Υποστύλωματα
Bc = 0.5  # Πλάτος (m)
Hc = 0.5  # Ύψος (m)

# Φορτία
hep = 0.05  # Πάχος επίστρωσης    (m)
htoix = 2.5 # Ύψος τοιχοπληρώσεων (m)

# Ειδικά βάρη 
gskyr = 25   # kN/m3
gepist = 20  # kN/m3
gtoix = 3.6  # kN/m2

# Φορτία Δοκού

# Μόνιμα φορτία πλακών

g1 = 1.5          # kN/m2
gIB = gskyr * hf  # kN/m2           
gol = g1 + gIB    # kN/m2

gd = gol * (a/2 + a/2)    # kN/m
gdtoix = gtoix * htoix    # kN/m
gdIB = gskyr * bw * (hb-hf)
g = gd + gdtoix + gdIB    # kN/m


# Ωφέλιμα φορτία πλακών
Q = 5                    # kN/m2
q = Q * (a/2 + a/2)      # kN/m


# Υλικά - Σκυρόδεμα C20/25 - Χάλυβας B500C

# Σκυρόδεμα 
fck = 20                            # Χαρακτηριστική Θλιπτική αντοχή (Mpa)
fcm = fck + 8                       # Μέση θλιπτική αντοχή (Mpa)
Ecm = round(22*(fcm/10)**0.3, 1)    # Μέτρο ελαστικότητας (Gpa) 
U = 0.0                             # Συντελεστής Poisson
E = Ecm * 10**6                     # GPa -> kN/m²
G = E / (2*(1+U))                   # Μέτρο διάτμησης (kN/m²)

# Διατομές

# Δοκός 
tbeam = TeeSectionGeometry(bw = bw, h = hb, beff=beff, hf = hf)
A_tbeam = tbeam.area
Iz_tbeam = tbeam.moment_of_inertia_xx * 0.5
Avy_tbeam = tbeam.shear_area_2 * 0.5

# Υποστύλωματα 
rect_col = RectangularSectionGeometry(b=Bc, h=Hc)
A_col = rect_col.area
Iz_col = rect_col.moment_of_inertia_xx * 0.5
Avy_col = rect_col.shear_area_2 * 0.5

# Μάζα
mass = (g+0.3*q)*L / 9.81

# Μοντελοποίηση 
ops.wipe()
ops.model('basic', '-ndm', 2, '-ndf', 3)

# Nodes 
ops.node(1, 0, 0)
ops.node(2, 0, H)
ops.node(3, L, 0)
ops.node(4, L, H)

def plotStructure(title):
    opsv.plot_model(fig_wi_he=(50,20))

    plt.title(title)
    plt.xlabel('Distance (m)')
    plt.ylabel('Distance (m)')
    plt.grid()
    plt.show()

plotStructure('Frame nodes')

# Στηρίξεις
ops.fix(1, 1, 1, 1)
ops.fix(3, 1, 1, 1)

ops.geomTransf('Linear', 1) 

ops.mass(int(2), mass, 1.0e-10, 1.0e-10)

ops.equalDOF(2, 4, 1)  

# Elements
# Υποστύλωματα
ops.element('ElasticTimoshenkoBeam', 1, 1, 2, E, G, A_col, Iz_col, Avy_col, 1)
ops.element('ElasticTimoshenkoBeam', 2, 3, 4, E, G, A_col, Iz_col, Avy_col, 1)
# Δοκός 
ops.element('ElasticTimoshenkoBeam', 3, 2, 4, E, G, A_tbeam, Iz_tbeam, Avy_tbeam, 1)

elem_type = {1:'Column', 2:'Column', 3:'Beam'}

plotStructure('Frame elements - Fixities')

numEigen = 1
eigenValues = ops.eigen('-genBandArpack', numEigen)

_periods = []
for i in range(0, numEigen):
    lamb = eigenValues[i]
    period = 2 * np.pi / np.sqrt(lamb)
    _periods.append(period)
    print(f'Period {i+1} = {period:.4f}s')

# Σεισμική δράση - Φάσμα σχεδιασμού EC8

agR = 0.36
ductility_class = 'M'
ground_type = 'C'
importance = 'II' 
γI=1.0
q_factor = 3.0 * 1.3
ag = agR*γI
specEC8 = SpectraEc8(αgR=agR,
                     γI=γI,
                     ground_type = ground_type,
                     spectrum_type = 1,
                     η=1.0,
                     q=q_factor,
                     β=0.2)

# Υπολογισμός σεισμικών δράσεων 

Sd_T = specEC8.Sd(period) * 9.81
M = mass 
if period <= 2*specEC8.TC:
    λ_ec8 = 0.85  # Για κτίρια με >2 ορόφους
else:
    λ_ec8 = 1.0   # Για μονοόροφα/μικρά

Vb = M * Sd_T * λ_ec8

print(f'Επιτάχυνση σχεδιασμού Sd(T) = {Sd_T:.3f}m/sec2')
print(f'Μάζα για το σύνολο του Φορέα M = {M:.2f}t')
print(f'Τέμνουσα βάσης Vb = {Vb:.2f}kN')

# Φορτία - Συνδυασμός 1: (G + 0.3Q + Ex)
ops.timeSeries('Constant', 1)
ops.pattern('Plain', 1, 1)

# Συγκεντρωμένο οριζόντιο φορτίο Vb στον κόμβο 2
ops.load(2, Vb, 0.0, 0.0)

# Ομοιόμορφο φορτίο μόνο στη δοκό (στοιχείο 3)
ops.eleLoad('-ele', 3, '-type', '-beamUniform', -(g + 0.3*q))

# Ανάλυση
ops.system('BandGeneral')
ops.numberer('RCM')
ops.constraints('Transformation')
ops.test('NormDispIncr', 1.0e-8, 6, 2)
ops.integrator('LoadControl', 1.0)
ops.algorithm('Linear')
ops.analysis('Static')

# Διαγράμματα 

# 1 step analysis
ops.analyze(1)
sfac = 200
opsv.plot_defo(sfac,
                   fig_wi_he=(50,20),
                   fmt_defo={'color': 'red', 'linestyle': (0, (4, 5)), 'linewidth': 1.5},
                   fmt_undefo={'color': 'green', 'linestyle': 'solid', 'linewidth': 2,},
                  )
plt.title('Deflection - G + 0.3Q + Ex Load Case')
plt.xlabel('Distance (m)')
plt.ylabel('Distance (m)')
plt.grid()
plt.show()

# Διαγράμματα Μ, V, N
sfacN, sfacV, sfacM = 7.e-3, 6.e-3, 4.e-3

opsv.section_force_diagram_2d('M', sfacM, fig_wi_he=(50,20),
                             fmt_secforce1={'color': 'green'},
                             fmt_secforce2={'color': 'green'})

plt.title('Bending Moment Diagram - G + 0.3Q + Ex')
plt.xlabel('Distance (m)')
plt.ylabel('Distance (m)')
plt.grid()
plt.show()

opsv.section_force_diagram_2d('V', sfacV, fig_wi_he=(50,20),
                             fmt_secforce1={'color': 'red'},
                             fmt_secforce2={'color': 'red'})

plt.title('Shear Force Diagram - G + 0.3Q + Ex')
plt.xlabel('Distance (m)')
plt.ylabel('Distance (m)')
plt.grid()
plt.show()

opsv.section_force_diagram_2d('N', sfacN, fig_wi_he=(50,20))

plt.title('Axial Force Diagram - G + 0.3Q + Ex')
plt.xlabel('Distance (m)')
plt.ylabel('Distance (m)')
plt.grid()
plt.show()

# Αφαίρεση προηγούμενου load pattern
ops.remove('loadPattern', 1)
ops.remove('timeSeries', 1)

# Μηδενισμός του domain (μετακινήσεις = 0)
ops.setTime(0.0)
ops.wipeAnalysis()

# Φορτία - Συνδυασμός 2: (1.35G + 1.5Q)
ops.timeSeries('Constant', 2)
ops.pattern('Plain', 2, 2)

ops.eleLoad('-ele', 3, '-type', '-beamUniform', -(1.35*g + 1.5*q))

# Ανάλυση 1.35G + 1.5Q
ops.system('BandGeneral')
ops.numberer('RCM')
ops.constraints('Transformation')
ops.test('NormDispIncr', 1.0e-8, 6, 2)
ops.integrator('LoadControl', 1.0)
ops.algorithm('Linear')
ops.analysis('Static')

# Διαγράμματα 
ops.analyze(1)
sfac = 100
opsv.plot_defo(sfac,
                   fig_wi_he=(50,20),
                   fmt_defo={'color': 'red', 'linestyle': (0, (4, 5)), 'linewidth': 1.5},
                   fmt_undefo={'color': 'green', 'linestyle': 'solid', 'linewidth': 2,},
                  )
plt.title('Deflection - 1.35G + 1.5Q Load Case')
plt.xlabel('Distance (m)')
plt.ylabel('Distance (m)')
plt.grid()
plt.show()


opsv.section_force_diagram_2d('M', sfacM, fig_wi_he=(50,20),
                             fmt_secforce1={'color': 'green'},
                             fmt_secforce2={'color': 'green'})

plt.title('Bending Moment Diagram - 1.35G + 1.5Q')
plt.xlabel('Distance (m)')
plt.ylabel('Distance (m)')
plt.grid()
plt.show()

opsv.section_force_diagram_2d('V', sfacV, fig_wi_he=(50,20),
                             fmt_secforce1={'color': 'red'},
                             fmt_secforce2={'color': 'red'})

plt.title('Shear Force Diagram - 1.35G + 1.5Q')
plt.xlabel('Distance (m)')
plt.ylabel('Distance (m)')
plt.grid()
plt.show()

opsv.section_force_diagram_2d('N', sfacN, fig_wi_he=(50,20))

plt.title('Axial Force Diagram - 1.35G + 1.5Q')
plt.xlabel('Distance (m)')
plt.ylabel('Distance (m)')
plt.grid()
plt.show()
