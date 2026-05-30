import openseespy.opensees as ops
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# ==========================================================
# 1. INITIALIZATION & FRAME GEOMETRICS (SI Units: N, mm)
# ==========================================================
ops.wipe()
ops.model('basic', '-ndm', 2, '-ndf', 3) 

H_story = 3000.0                 # Story Height
L_bay = 6000.0                   # Bay Width
B_col, D_col = 400.0, 400.0      # Column Dimensions
B_beam, D_beam = 400.0, 600.0    # Beam Dimensions
Ag_col = B_col * D_col   

# Node Grid Definition
ops.node(1, 0.0,   0.0)
ops.node(2, L_bay, 0.0)
ops.node(3, 0.0,   H_story)
ops.node(4, L_bay, H_story)
ops.node(5, 0.0,   2.0 * H_story)
ops.node(6, L_bay, 2.0 * H_story)

# Fixed Boundary Conditions
ops.fix(1, 1, 1, 1)
ops.fix(2, 1, 1, 1)

# ==========================================================
# 2. MATERIAL CONSTITUTIVE LAWS (MPa)
# ==========================================================
# Concrete Core & Cover Parameters
fcc, ecc0, fccu, eccu = -48.0, -0.004, -10.0, -0.0125       
fc, ec0, fcu, ecu = -40.0, -0.002, 0.0, -0.005
ft, Ets = 4.0, 4000.0             

ops.uniaxialMaterial('Concrete02', 1, fcc, ecc0, fccu, eccu, 0.1, ft, Ets) # Confined Core
ops.uniaxialMaterial('Concrete02', 2, fc, ec0, fcu, ecu, 0.1, ft, Ets)     # Unconfined Cover

# Steel Reinforcement Parameters (Explicit Variable Binding)
fy = 400.0                       # Yield Strength (MPa) - Fixed NameError
Es = 200000.0                    # Elastic Modulus (MPa)
b_ratio = 0.015                  # Strain Hardening Ratio
R0, cR1, cR2 = 18.0, 0.925, 0.15 

ops.uniaxialMaterial('Steel02', 3, fy, Es, b_ratio, R0, cR1, cR2)          # Rebar Steel

# ==========================================================
# 3. FIBER SECTION GENERATION (Column vs Beam Discretization)
# ==========================================================
cover = 40.0        
bar_area = np.pi * (20.0**2) / 4.0  

# --- Section 1: Column (400x400 mm) ---
core_Y_col = (B_col / 2.0) - cover   
core_Z_col = (D_col / 2.0) - cover   
ops.section('Fiber', 1)
ops.patch('quad', 1, 16, 16, -core_Y_col, -core_Z_col, core_Y_col, -core_Z_col, core_Y_col, core_Z_col, -core_Y_col, core_Z_col)
ops.patch('quad', 2, 4, 16, -B_col/2, -D_col/2,  core_Y_col, -D_col/2,  core_Y_col, -core_Z_col, -B_col/2, -core_Z_col)
ops.patch('quad', 2, 4, 16, -core_Y_col,  core_Z_col,   B_col/2,  core_Z_col,   B_col/2,  D_col/2, -core_Y_col,  D_col/2)
ops.patch('quad', 2, 16, 4, -B_col/2, -core_Z_col, -core_Y_col, -core_Z_col, -core_Y_col,  core_Z_col, -B_col/2,  core_Z_col)
ops.patch('quad', 2, 16, 4,  core_Y_col,  -core_Z_col,  B_col/2, -core_Z_col,  B_col/2,  core_Z_col,  core_Y_col,   core_Z_col)
for by, bz in [(-core_Y_col, -core_Z_col), (0, -core_Z_col), (core_Y_col, -core_Z_col), (-core_Y_col, 0), (core_Y_col, 0), (-core_Y_col, core_Z_col), (0, core_Z_col), (core_Y_col, core_Z_col)]:
    ops.fiber(by, bz, bar_area, 3)

# --- Section 2: Beam (400x600 mm) ---
core_Y_beam = (D_beam / 2.0) - cover   
core_Z_beam = (B_beam / 2.0) - cover   
ops.section('Fiber', 2)
ops.patch('quad', 1, 16, 16, -core_Y_beam, -core_Z_beam, core_Y_beam, -core_Z_beam, core_Y_beam, core_Z_beam, -core_Y_beam, core_Z_beam)
ops.patch('quad', 2, 4, 16, -D_beam/2, -B_beam/2,  core_Y_beam, -B_beam/2,  core_Y_beam, -core_Z_beam, -D_beam/2, -core_Z_beam)
ops.patch('quad', 2, 4, 16, -core_Y_beam,  core_Z_beam,   D_beam/2,  core_Z_beam,   D_beam/2,  B_beam/2, -core_Y_beam,  B_beam/2)
ops.patch('quad', 2, 16, 4, -D_beam/2, -core_Z_beam, -core_Y_beam, -core_Z_beam, -core_Y_beam,  core_Z_beam, -D_beam/2,  core_Z_beam)
ops.patch('quad', 2, 16, 4,  core_Y_beam,  -core_Z_beam,  D_beam/2, -core_Z_beam,  D_beam/2,  core_Z_beam,  core_Y_beam,   core_Z_beam)
for by, bz in [(-core_Y_beam, -core_Z_beam), (0, -core_Z_beam), (core_Y_beam, -core_Z_beam), (-core_Y_beam, 0), (core_Y_beam, 0), (-core_Y_beam, core_Z_beam), (0, core_Z_beam), (core_Y_beam, core_Z_beam)]:
    ops.fiber(by, bz, bar_area, 3)

ops.geomTransf('Linear', 1)

# Element Connectivity Matrix
ops.element('nonlinearBeamColumn', 1, 1, 3, 5, 1, 1) # Lower Left Col
ops.element('nonlinearBeamColumn', 2, 2, 4, 5, 1, 1) # Lower Right Col
ops.element('nonlinearBeamColumn', 3, 3, 5, 5, 1, 1) # Upper Left Col
ops.element('nonlinearBeamColumn', 4, 4, 6, 5, 1, 1) # Upper Right Col
ops.element('nonlinearBeamColumn', 5, 3, 4, 5, 2, 1) # Floor Beam
ops.element('nonlinearBeamColumn', 6, 5, 6, 5, 2, 1) # Roof Beam

# ==========================================================
# 4. GRAVITY LOADING & ANALYSIS SETUP
# ==========================================================
P_axial = -0.10 * 40.0 * Ag_col  
ops.timeSeries('Linear', 1)
ops.pattern('Plain', 1, 1)
ops.load(5, 0.0, P_axial, 0.0)
ops.load(6, 0.0, P_axial, 0.0)

ops.constraints('Transformation')
ops.numberer('RCM')
ops.system('BandGeneral')
ops.test('NormDispIncr', 1.0e-6, 30)
ops.algorithm('Newton')
ops.integrator('LoadControl', 0.1)
ops.analysis('Static')
ops.analyze(10)
ops.loadConst('-time', 0.0)

# ==========================================================
# 5. CYCLIC RUN CONTROL (NODE 5 ROOF CONTROL)
# ==========================================================
drift_protocols = [0.005, 0.010, 0.015, 0.020, 0.030] 
cycles_per_level = 3
step_size = 0.00025                                   
H_total = 2.0 * H_story

ops.pattern('Plain', 2, 1)
ops.load(5, 1.0, 0.0, 0.0)

all_drifts, all_shears = [], []

for drift in drift_protocols:
    peak_disp = drift * H_total 
    for cycle in range(cycles_per_level):
        ops.integrator('DisplacementControl', 5, 1, step_size * H_total)
        ops.analysis('Static')
        while ops.nodeDisp(5, 1) < peak_disp:
            if ops.analyze(1) != 0: break
            all_drifts.append(ops.nodeDisp(5, 1) / H_total * 100)
            ops.reactions()
            all_shears.append(-(ops.nodeReaction(1, 1) + ops.nodeReaction(2, 1)) / 1000) 
            
        ops.integrator('DisplacementControl', 5, 1, -step_size * H_total)
        while ops.nodeDisp(5, 1) > -peak_disp:
            if ops.analyze(1) != 0: break
            all_drifts.append(ops.nodeDisp(5, 1) / H_total * 100)
            ops.reactions()
            all_shears.append(-(ops.nodeReaction(1, 1) + ops.nodeReaction(2, 1)) / 1000)
            
        ops.integrator('DisplacementControl', 5, 1, step_size * H_total)
        while ops.nodeDisp(5, 1) < 0.0:
            if ops.analyze(1) != 0: break
            all_drifts.append(ops.nodeDisp(5, 1) / H_total * 100)
            ops.reactions()
            all_shears.append(-(ops.nodeReaction(1, 1) + ops.nodeReaction(2, 1)) / 1000)

# ==========================================================
# 6. PLASTIC CAPACITY VALIDATION METRICS
# ==========================================================
Ast_tension = 3.0 * bar_area  
Mp_col = ((Ast_tension * fy * (2.0 * core_Y_col)) + (abs(P_axial) * (B_col / 2.0 - cover / 2.0))) / 1.0e6  
Mp_beam = (Ast_tension * fy * (2.0 * core_Y_beam)) / 1.0e6  
Vn_frame_theoretical = (Mp_col + 2.0 * Mp_beam) / (H_story / 1000.0)

# ==========================================================
# 7. VISUALIZATION PLATFORM WITH DIMENSION METRICS
# ==========================================================
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 9.5
fig = plt.figure(figsize=(16, 9.2))
gs = fig.add_gridspec(2, 3, width_ratios=[1.2, 1.0, 1.8], height_ratios=[4.5, 1.0])

# --- Subplot 1: Dimensioned Kinematic Model Layout ---
ax_schem = fig.add_subplot(gs[0, 0])
ax_schem.plot([0, 0], [0, 2*H_story], color='#4b5563', linewidth=10, solid_capstyle='butt', zorder=1)
ax_schem.plot([L_bay, L_bay], [0, 2*H_story], color='#4b5563', linewidth=10, solid_capstyle='butt', zorder=1)
ax_schem.plot([0, L_bay], [H_story, H_story], color='#1e3a8a', linewidth=14, solid_capstyle='butt', zorder=1) 
ax_schem.plot([0, L_bay], [2*H_story, 2*H_story], color='#1e3a8a', linewidth=14, solid_capstyle='butt', zorder=1) 

# Boundaries & Footings
ax_schem.plot([-500, 500], [0, 0], color='#111111', linewidth=4, zorder=2)
ax_schem.plot([L_bay-500, L_bay+500], [0, 0], color='#111111', linewidth=4, zorder=2)
for bx in [0, L_bay]:
    for k in np.linspace(bx-400, bx+400, 6):
        ax_schem.plot([k, k-100], [0, -150], color='#111111', linewidth=1.5)

# Clear System Dimension Lines
ax_schem.annotate('', xy=(0, -350), xytext=(L_bay, -350), arrowprops=dict(arrowstyle='<->', color='#334155', lw=1.2))
ax_schem.text(L_bay/2, -500, f'L = {L_bay:.0f} mm', ha='center', va='top', color='#334155', weight='bold', fontsize=8.5)

ax_schem.annotate('', xy=(-600, 0), xytext=(-600, H_story), arrowprops=dict(arrowstyle='<->', color='#334155', lw=1.2))
ax_schem.text(-750, H_story/2, f'H1 = {H_story:.0f} mm', ha='right', va='center', color='#334155', weight='bold', fontsize=8.5)

ax_schem.annotate('', xy=(-600, H_story), xytext=(-600, 2*H_story), arrowprops=dict(arrowstyle='<->', color='#334155', lw=1.2))
ax_schem.text(-750, 1.5*H_story, f'H2 = {H_story:.0f} mm', ha='right', va='center', color='#334155', weight='bold', fontsize=8.5)

# Loading Vector Labels
ax_schem.annotate('', xy=(0, 2*H_story), xytext=(-1200, 2*H_story), arrowprops=dict(facecolor='#1e40af', edgecolor='#1e40af', width=2.5, headwidth=8))
ax_schem.text(-1250, 2*H_story, 'Cyclic Actuator\nInput (3% Drift)', color='#1e40af', weight='bold', fontsize=8.5, ha='right', va='center')

ax_schem.set_xlim([-2400, L_bay + 1200])
ax_schem.set_ylim([-700, 2*H_story + 1200])
ax_schem.set_title("2-Story Frame Kinematics & Dimensions", fontsize=10.5, fontweight='bold', pad=12)
ax_schem.axis('off')

# --- Subplot 2: Multi-Section Profile Layout with Dimensions ---
gs_sec = gs[0, 1].subgridspec(2, 1, hspace=0.38)
ax_col_sec = fig.add_subplot(gs_sec[0, 0])
ax_beam_sec = fig.add_subplot(gs_sec[1, 0])

# Column Section Render
ax_col_sec.add_patch(patches.Rectangle((-200, -200), 400, 400, facecolor='#f8fafc', edgecolor='#0f172a', linewidth=1.5))
ax_col_sec.add_patch(patches.Rectangle((-160, -160), 320, 320, facecolor='#e2e8f0', edgecolor='#b91c1c', linestyle='--', linewidth=1.2))
for by, bz in [(-160,-160), (0,-160), (160,-160), (-160,0), (160,0), (-160,160), (0,160), (160,160)]:
    ax_col_sec.add_patch(patches.Circle((by, bz), 14, facecolor='#0f172a', zorder=4))
ax_col_sec.annotate('', xy=(-200, -235), xytext=(200, -235), arrowprops=dict(arrowstyle='<->', color='#475569', lw=1))
ax_col_sec.text(0, -255, "400 mm", ha='center', va='top', fontsize=8, color='#475569', weight='bold')
ax_col_sec.set_xlim([-260, 260]); ax_col_sec.set_ylim([-260, 260]); ax_col_sec.set_aspect('equal')
ax_col_sec.set_title("Column Profile (400x400 mm)", fontsize=9, fontweight='bold')
ax_col_sec.grid(True, linestyle=':', alpha=0.4)

# Beam Section Render
ax_beam_sec.add_patch(patches.Rectangle((-300, -200), 600, 400, facecolor='#f8fafc', edgecolor='#0f172a', linewidth=1.5))
ax_beam_sec.add_patch(patches.Rectangle((-260, -160), 520, 320, facecolor='#cbd5e1', edgecolor='#b91c1c', linestyle='--', linewidth=1.2))
for by, bz in [(-260,-160), (0,-160), (260,-160), (-260,0), (260,0), (-260,160), (0,160), (260,160)]:
    ax_beam_sec.add_patch(patches.Circle((by, bz), 14, facecolor='#0f172a', zorder=4))
ax_beam_sec.annotate('', xy=(-300, -240), xytext=(300, -240), arrowprops=dict(arrowstyle='<->', color='#475569', lw=1))
ax_beam_sec.text(0, -260, "600 mm (Depth)", ha='center', va='top', fontsize=8, color='#475569', weight='bold')
ax_beam_sec.annotate('', xy=(-340, -200), xytext=(-340, 200), arrowprops=dict(arrowstyle='<->', color='#475569', lw=1))
ax_beam_sec.text(-360, 0, "400 mm", ha='right', va='center', fontsize=8, color='#475569', weight='bold', rotation=90)
ax_beam_sec.set_xlim([-390, 390]); ax_beam_sec.set_ylim([-270, 270]); ax_beam_sec.set_aspect('equal')
ax_beam_sec.set_title("Beam Profile (400x600 mm)", fontsize=9, fontweight='bold')
ax_beam_sec.grid(True, linestyle=':', alpha=0.4)

# --- Subplot 3: Global Response Curve ---
ax_curve = fig.add_subplot(gs[0, 2])
ax_curve.plot(all_drifts, all_shears, color='#1e40af', linewidth=1.5, alpha=0.95, label='OpenSees Frame System Response')
ax_curve.axhline(Vn_frame_theoretical, color='#dc2626', linestyle='--', linewidth=1.8, label=f'Theoretical Capacity $V_n$ = $\pm${Vn_frame_theoretical:.1f} kN')
ax_curve.axhline(-Vn_frame_theoretical, color='#dc2626', linestyle='--', linewidth=1.8)
ax_curve.grid(True, linestyle='--', alpha=0.5)
ax_curve.set_xlim([-3.5, 3.5]); ax_curve.set_ylim([-350, 350])
ax_curve.set_title("System Hysteresis Performance Curve", fontsize=10.5, fontweight='bold', pad=12)
ax_curve.set_xlabel("Global Roof Drift Ratio (%)", fontweight='bold', fontsize=9)
ax_curve.set_ylabel("Total Base Shear Force $V_b$ (kN)", fontweight='bold', fontsize=9)
ax_curve.legend(loc='upper left', fontsize=8.5, framealpha=1.0)

# --- Row 2: Comprehensive Structural Calculation Log Block ---
ax_summary = fig.add_subplot(gs[1, :])
ax_summary.axis('off')  
summary_text = (
    "=================================== REINFORCED CONCRETE MULTI-MEMBER FRAME REPORT ===================================\n"
    f"  • Frame Topology: 2-Story 1-Bay Frame Layout  |  Bay Width = {L_bay:.0f} mm  |  Story Height = {H_story:.0f} mm (Total Height = {2*H_story:.0f} mm)\n"
    f"  • Member Configurations: Columns = 400x400 mm (8d20)  |  Beams = 400x600 mm (8d20, Extended Lever Arm d' = {2*core_Y_beam:.0f} mm)\n"
    f"  • Yield Mechanics (First-Principles Ultimate Evaluation):\n"
    f"      M_p_column (at P = 640 kN) = {Mp_col:.1f} kN·m   |   M_p_beam (Pure Bending) = {Mp_beam:.1f} kN·m\n"
    f"      Governing Strong-Column/Weak-Beam Dissipation Mechanism Code Calculation Formulation:\n"
    f"      V_n = (M_p_col + 2 * M_p_beam) / H_story = ({Mp_col:.1f} + 2 * {Mp_beam:.1f}) / {H_story/1000.0:.1f} = {Vn_frame_theoretical:.1f} kN"
)
ax_summary.text(0.5, 0.4, summary_text, fontsize=8.5, family='monospace', weight='bold', color='#0f172a',
                bbox=dict(boxstyle='round,pad=0.5', facecolor='#f8fafc', edgecolor='#1e40af', linewidth=1.5), ha='center', va='center')

plt.suptitle("OpenSeesPy Seismic Behavior Evaluation Dashboard for Multi-Story RC Frames", fontsize=13, fontweight='bold', y=0.97)
plt.subplots_adjust(left=0.05, right=0.96, top=0.88, bottom=0.05, wspace=0.28)
plt.show()