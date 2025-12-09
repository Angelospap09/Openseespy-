# Openseespy Timoshenko Beam Element Example

import openseespy.opensees as ops
import numpy as np
import matplotlib.pyplot as plt
import opsvis as opsv
import math 
from streng.ppp.sections.geometry.tee import TeeSectionGeometry
from streng.ppp.sections.geometry.rectangular import RectangularSectionGeometry
from streng.codes.eurocodes.ec2.raw.ch5.geometric_data import effective_width
from streng.codes.eurocodes.ec8.cls.seismic_action.spectra import SpectraEc8

# Δεδομένα

# Μονάδες: kN, m, C

# Διαστάσεις

L = 7             # Μήκος δοκού (m)
H = 3             # Mήκος υποστυλώματος (m)
a = L/4           # Απόσταση φορτίου από ακραίο υποστύλωμα (m)
L1 = L/2          # Μήκος τμήματος δοκού μεταξύ υποστυλωμάτων (m)
cnom = 0.05       # Επικάλυψη (m)

# Γεωμετρικοί παραμέτροι 
nBays = 2       # Αριθμός ανοιγμάτων 
nStories = 3    # Αριθμός ορόφων
h_bay = H       # Ύψος ορόφου (m)
w_bay = L       # Πλάτος ανοίγματος (m)

# Πάκτωση
fixity = (1,1,1)

# Δοκός - Τυπική διατομή Τ
bw = 0.3         
hf = 0.15        
hb = 0.5         

# zero_moments_case: 1=ακραίο άνοιγμα (l0=0.85*L), 2=μεσαίο (l0=0.7*L), 3=πρόβολος
zero_moments_case = 1 if nBays == 1 else 2
l0 = effective_width.l0(l1=L, l2=L, zero_moments_case=zero_moments_case)
b1 = b2 = L/2 - bw/2   
_b = bw + b1 + b2
beff1 = beff2 = effective_width.beffi(b1, l0)
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
print(f'g = {g:.3f} kN/m')

# Ωφέλιμα φορτία πλακών
Q = 5                    # kN/m2
q = Q * (a/2 + a/2)      # kN/m
print(f'q = {q:.3f} kN/m')

# Υλικά - Σκυρόδεμα C20/25 - Χάλυβας B500C

# Σκυρόδεμα 
fck = 20                            # Χαρακτηριστική Θλιπτική αντοχή (Mpa)
fcm = fck + 8                       # Μέση θλιπτική αντοχή (Mpa)
Ecm = round(22*(fcm/10)**0.3, 1)    # Μέτρο ελαστικότητας (GPa) - EC2 τύπος
U = 0.0                             # Συντελεστής Poisson
E = Ecm * 1e6                       # GPa -> kN/m² 
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

# Μάζα (ανά όροφο, λαμβάνοντας όλα τα ανοίγματα)
mass = nBays * (g + 0.3*q) * L / 9.81  # kN → t για κάθε όροφο

# Συγκέντρωση μαζών ανά όροφο (master node κάθε ορόφου)
story_masses = []
for story in range(nStories):
    story_masses.append(mass)

# Model Opensees
ops.wipe()
ops.model('basic', '-ndm', 2, '-ndf', 3)


# Nodes - Σωστή αρίθμηση: κάθε όροφο έχει 2 κόμβους ανά στήλη (κάτω και πάνω)
for story in range(nStories + 1):  
    for i in range(nBays + 1):     
        index = i + 1
        xCoord = i * w_bay
        yCoord = story * h_bay
        
        # Κόμβος σε αυτό το ύψος
        nodeTag = story * (nBays + 1) + i + 1
        ops.node(nodeTag, xCoord, yCoord)
        print(f'Node {nodeTag} at ({xCoord:.1f}, {yCoord:.1f}) defined')




def plotStructure(title):
    #Pass in optional argument (fig_wi_he) to control plot size
    opsv.plot_model(fig_wi_he=(50,20))

    plt.title(title)
    plt.xlabel('Distance (m)')
    plt.ylabel('Distance (m)')
    plt.grid()
    plt.show()

plotStructure('Frame nodes')

transftype = 'Linear'
transftag = 1
ops.geomTransf(transftype, transftag)

# Define Υποστυλώματα 
eletag = 1
for story in range(nStories):
    for i in range(nBays + 1):
        nodeBottom = story * (nBays + 1) + i + 1
        nodeTop = (story + 1) * (nBays + 1) + i + 1
        
        ops.element('ElasticTimoshenkoBeam', eletag, nodeBottom, nodeTop, E, G, A_col, Iz_col, Avy_col, transftag)
        print(f'Column element {eletag}: nodes {nodeBottom}-{nodeTop}')
        eletag += 1

# Define Δοκοί - Loop για κάθε όροφο
eletag_beam = eletag  # Συνέχιση από τα columns
for story in range(nStories):
    for i in range(nBays):
        # Nodes στην κορυφή του ορόφου (story+1)
        nodeLeft = (story + 1) * (nBays + 1) + i + 1      # Αριστερή: story+1, col i
        nodeRight = (story + 1) * (nBays + 1) + (i + 1) + 1   # Δεξιά: story+1, col i+1
        
        ops.element('ElasticTimoshenkoBeam', eletag_beam, nodeLeft, nodeRight, E, G, A_tbeam, Iz_tbeam, Avy_tbeam, transftag)
        print(f'Beam element {eletag_beam}: nodes {nodeLeft}-{nodeRight}')
        eletag_beam += 1

plotStructure('Frame elements')

# Στηρίξεις - Μόνο η βάση (story 0)
for i in range(nBays + 1):
    nodeTag = 0 * (nBays + 1) + i + 1
    ops.fix(nodeTag, 1, 1, 1)

plotStructure('Frame supports')

# Ιδιοπερίοδος Κατασκευής

# Μάζα στον αριστερό κόμβο κάθε ορόφου (ΤΟΠ του ορόφου)
for story in range(nStories):
    nodeTag = (story + 1) * (nBays + 1) + 1  # Αριστερός κόμβος της κορυφής (story+1)
    ops.mass(nodeTag, mass, 1.0e-10, 1.0e-10)

# equalDOF για κάθε όροφο (κοινή οριζόντια μετακίνηση)
for story in range(1, nStories + 1):
    master = story * (nBays + 1) + 1
    for i in range(1, nBays + 1):
        slave = story * (nBays + 1) + i + 1
        ops.equalDOF(master, slave, 1)  # DOF 1 = X direction

numEigen = 1
eigenValues = ops.eigen('-genBandArpack', numEigen)

_periods = []
for i in range(0, numEigen):
    lamb = eigenValues[i]
    period = 2 * np.pi / np.sqrt(lamb)
    _periods.append(period)
    print(f'Period {i+1} = {period:.4f}s')

# Ιδιομορφικό διάνυσμα για κάθε master node (όροφοι 1, 2, 3)
eigen_vector = []
master_nodes = [(story + 1) * (nBays + 1) + 1 for story in range(nStories)]
for mn in master_nodes:
    eigen_vector.append(ops.nodeEigenvector(mn, 1, 1))

print('Συνιστώσες 1ης ιδιομορφής (δε θα δείτε ακριβώς τα ίδια νούμερα, θα έχετε όμως τις ίδιες αναλογίες):')
print([f'{e:.4f}' for e in eigen_vector])

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
M = mass * nStories  # συνολική μάζα φορέα

if period <= 2*specEC8.TC:
    λ_ec8 = 0.85
else:
    λ_ec8 = 1.0
    
Fb = M * Sd_T * λ_ec8

eigen_vector_array = np.array(eigen_vector)
mass_array = np.array([mass, mass, mass])

print(f'Επιτάχυνση σχεδιασμού Sd(T) = {Sd_T:.3f}m/sec2')
print(f'Μάζα για το σύνολο του κτιρίου M = {M:.2f}t')
print(f'λ = {λ_ec8}')
print(f'Τέμνουσα βάσης Fb = {Fb:.2f}kN')

# Φορτία - Συνδυασμός 1: (1.35G + 1.5Q)
ops.timeSeries('Constant', 1)
ops.pattern('Plain', 1, 1)

# Ομοιόμορφο φορτίο σε όλες τις δοκούς
eletag_beam = (nStories * (nBays + 1)) + 1  # Πρώτο beam element
for story in range(nStories):
    for i in range(nBays):
        ops.eleLoad('-ele', eletag_beam, '-type', '-beamUniform', -(1.35*g + 1.5*q))
        eletag_beam += 1

# Διάγραμμα φορτίων με χειροκίνητη προσθήκη κειμένων
def plotLoads(title, show_seismic=True, dist_value=None, dist_label=None, Fi_plot=None):
    fig, ax = plt.subplots(figsize=(12, 8))
    
    try:
        # Σχεδίαση φορτίων/στηρίξεων από το ενεργό load pattern
        opsv.plot_load(
            nep=11,
            sfac=False,
            fig_wi_he=False,
            fig_lbrt=False,
            fmt_model_loads={'color': 'red', 'linestyle': 'solid', 'linewidth': 1.5, 'marker': '', 'markersize': 1},
            node_supports=True,
            truss_node_offset=0,
            local_axes=False,
            ax=ax
        )

        # Προσθήκη annotations για οριζόντια Fi (αν δοθούν)
        if show_seismic and Fi_plot is not None:
            for story, fi_val in enumerate(Fi_plot):
                y_pos = (story + 1) * h_bay
                ax.annotate(f'Fi = {fi_val:.1f} kN',
                           xy=(-0.5, y_pos),
                           fontsize=10,
                           color='blue',
                           weight='bold',
                           bbox=dict(boxstyle='round,pad=0.3', facecolor='lightblue', alpha=0.7))
                ax.arrow(-1, y_pos, 0.8, 0, head_width=0.2, head_length=0.15,
                         fc='blue', ec='blue', linewidth=2)

        # Ετικέτες για κατανεμημένα
        dist_val = dist_value if dist_value is not None else (g + 0.3*q)
        dist_lab = dist_label if dist_label is not None else 'g + 0.3q'
        for story in range(nStories):
            y_pos = (story + 1) * h_bay
            x_mid = w_bay / 2
            ax.text(x_mid, y_pos + 0.3, f'{dist_lab} = {dist_val:.1f} kN/m',
                    fontsize=9, ha='center', color='red',
                    bbox=dict(boxstyle='round,pad=0.3', facecolor='pink', alpha=0.7))

        ax.set_title(title, fontsize=14, weight='bold')
        ax.set_xlabel('Distance (m)', fontsize=11)
        ax.set_ylabel('Distance (m)', fontsize=11)
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.show()
    except Exception as e:
        print(f"Could not plot loads: {e}")

plotLoads('Loading: 1.35G + 1.5Q', show_seismic=False, dist_value=(1.35*g + 1.5*q), dist_label='1.35G + 1.5Q', Fi_plot=None)

# Ανάλυση

# Create SOE
ops.system('BandGeneral')

# Create DOF number
ops.numberer('RCM')

# Create constraint handler
ops.constraints('Transformation')

# Create integrator
ops.integrator('LoadControl', 1)

# Create convergence test
ops.test('NormDispIncr', 1.0e-8, 6, 2)

# Create algorithm
ops.algorithm('Linear')

# Create analysis object
ops.analysis('Static')

# 1 step analysis
ops.analyze(1)
sfac = 200

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

# Διαγράμματα Μ, V, N
sfacN, sfacV, sfacM = 7.e-3, 6.e-3, 4.e-3

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

# Μετακινήσεις
u2 = ops.nodeDisp(2, 2)
u4 = ops.nodeDisp(4, 2)

# Αφαίρεση προηγούμενου load pattern
ops.remove('loadPattern', 1)
ops.remove('timeSeries', 1)
ops.wipeAnalysis()
ops.loadConst('-time', 0.0)
ops.setTime(0.0)
ops.wipeAnalysis()

# Φορτία - Συνδυασμός 2: (G + 0.3Q + Ex) με κατανομή Fi
ops.timeSeries('Constant', 2)
ops.pattern('Plain', 2, 2)

# Κατανομή σεισμικών δυνάμεων ανά όροφο με βάση την ιδιομορφή
Fi = (eigen_vector_array * mass_array) * Fb / sum(eigen_vector_array * mass_array)
for _i, _fi in enumerate(Fi):
    print(f'F{_i+1} = {_fi:.2f}kN')

# Συγκεντρωμένα οριζόντια Fi στα master nodes κάθε ορόφου
for story, fi in enumerate(Fi):
    nodeTag = (story + 1) * (nBays + 1) + 1
    ops.load(nodeTag, fi, 0.0, 0.0)

# Ομοιόμορφο φορτίο σε όλες τις δοκούς
eletag_beam = (nStories * (nBays + 1)) + 1  # Πρώτο beam element
for story in range(nStories):
    for i in range(nBays):
        ops.eleLoad('-ele', eletag_beam, '-type', '-beamUniform', -(g + 0.3*q))
        eletag_beam += 1

# Διάγραμμα φορτίων για G + 0.3Q + Ex (με Fi)
plotLoads('Loading: G + 0.3Q + Ex', show_seismic=True, dist_value=(g + 0.3*q), dist_label='G + 0.3Q', Fi_plot=Fi)


# Ανάλυση G + 0.3Q + Ex
ops.system('BandGeneral')
ops.numberer('RCM')
ops.constraints('Transformation')
ops.test('NormDispIncr', 1.0e-8, 6, 2)
ops.integrator('LoadControl', 1.0)
ops.algorithm('Linear')
ops.analysis('Static')

# Διαγράμματα
ops.analyze(1)
sfac = 150
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



