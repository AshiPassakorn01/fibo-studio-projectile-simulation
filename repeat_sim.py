import tkinter as tk
from tkinter import ttk, messagebox
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import math
import tools
from physic import PhysicsEngine
import csv
import os
from datetime import datetime

class RepeatabilityWindow:
    def __init__(self):
        self.root = tk.Tk()
        self.root.title("Repeatability & Error Analysis (6 Cases + Analytical Bounds)")
        self.root.geometry("1650x900") 
        
        self.last_run_data = None 
        
        try:
            self.sd = tools.load_data()
        except:
            messagebox.showerror("Error", "Data file not found. Please run the Simulation from the main window first.")
            self.root.destroy()
            return
            
        self.setup_ui()
        
    def setup_ui(self):
        # ==========================================
        # 1. LEFT PANEL (Inputs & Base Params)
        # ==========================================
        left_frame = tk.Frame(self.root, width=340, bg="#f4f4f4", relief=tk.RIDGE, bd=2)
        left_frame.pack(side=tk.LEFT, fill=tk.Y, padx=5, pady=10)
        left_frame.pack_propagate(False)
        
        tk.Label(left_frame, text="⚙️ Error Parameters", font=("Arial", 12, "bold"), bg="#f4f4f4").pack(pady=10)
        
        self.entries = {}
        def add_input(label, key, default, unit):
            f = tk.Frame(left_frame, bg="#f4f4f4")
            f.pack(fill=tk.X, padx=10, pady=5)
            tk.Label(f, text=label, width=17, anchor="w", bg="#f4f4f4").pack(side=tk.LEFT)
            ent = tk.Entry(f, width=8, justify="center")
            ent.insert(0, default)
            ent.pack(side=tk.LEFT, padx=5)
            tk.Label(f, text=unit, bg="#f4f4f4").pack(side=tk.LEFT)
            self.entries[key] = ent

        add_input("Number of Shots (N):", "n_shots", "200", "shots")
        add_input("Target Tolerance:", "tolerance", "60", "mm")
        ttk.Separator(left_frame, orient='horizontal').pack(fill=tk.X, pady=10, padx=10)
        
        # แยก Pitch และ Yaw
        add_input("1. Pitch Error (±):", "err_pitch", "0.5", "deg")
        add_input("2. Yaw Error (±):", "err_yaw", "0.5", "deg")
        add_input("3. Pressure Error (±):", "err_press", "0.02", "bar")
        add_input("4. Delay Error (±):", "err_delay", "0.03", "sec")
        
        btn_frame = tk.Frame(left_frame, bg="#f4f4f4")
        btn_frame.pack(fill=tk.X, padx=15, pady=20)
        
        tk.Button(btn_frame, text="▶ Run Simulation", font=("Arial", 11, "bold"), bg="#0275d8", fg="white", height=2, command=self.run_simulation).pack(fill=tk.X, pady=(0, 10))
        self.btn_export = tk.Button(btn_frame, text="💾 Export Results", font=("Arial", 11, "bold"), bg="#5cb85c", fg="white", height=2, command=self.export_results, state=tk.DISABLED)
        self.btn_export.pack(fill=tk.X)

        base_frame = tk.LabelFrame(left_frame, text=" Base Parameters ", bg="#f4f4f4", font=("Arial", 10, "bold"))
        base_frame.pack(fill=tk.X, padx=10, pady=10)
        tk.Label(base_frame, text=f"Pitch: {self.sd.get('pitch',0):.2f}° | Yaw: {self.sd.get('yaw',0):.2f}°\n"
                                  f"Pressure: {self.sd.get('pressure',0)} bar\n"
                                  f"Initial Vel (V0): {self.sd.get('v0',0):.3f} m/s\n"
                                  f"Delay: {self.sd.get('t_delay',0):.3f} s", justify="left", bg="#f4f4f4").pack(padx=5, pady=10)

        # ==========================================
        # 2. RIGHT PANEL (Math Bounds & Results Table)
        # ==========================================
        right_frame = tk.Frame(self.root, width=420, bg="#ffffff", relief=tk.RIDGE, bd=2)
        right_frame.pack(side=tk.RIGHT, fill=tk.Y, padx=5, pady=10)
        right_frame.pack_propagate(False)

        # --- ส่วนแสดงสมการทางคณิตศาสตร์ ---
        self.math_frame = tk.LabelFrame(right_frame, text=" 📐 Analytical Boundaries ", bg="#ffffff", font=("Arial", 10, "bold"), fg="#8e44ad")
        self.math_frame.pack(fill=tk.X, padx=10, pady=10)
        self.lbl_math = tk.Label(self.math_frame, text="Press 'Run Simulation' to\ncalculate dynamic boundaries.", justify="left", bg="#ffffff", font=("Consolas", 9))
        self.lbl_math.pack(padx=5, pady=5, anchor="w")

        tk.Label(right_frame, text="📊 Statistical Results", font=("Arial", 12, "bold"), bg="#ffffff").pack(pady=(10, 5))
        
        cols = ("Case", "Rate", "MaxErr", "MeanErr")
        self.tree_res = ttk.Treeview(right_frame, columns=cols, show="headings", height=15)
        self.tree_res.heading("Case", text="Case")
        self.tree_res.heading("Rate", text="Hit Rate")
        self.tree_res.heading("MaxErr", text="Max Err(mm)")
        self.tree_res.heading("MeanErr", text="Mean Err(mm)")
        
        self.tree_res.column("Case", width=180, anchor="w")
        self.tree_res.column("Rate", width=70, anchor="center")
        self.tree_res.column("MaxErr", width=80, anchor="center")
        self.tree_res.column("MeanErr", width=80, anchor="center")
        self.tree_res.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)

        # ==========================================
        # 3. CENTER PANEL (Plots - 2x3 Grid)
        # ==========================================
        center_frame = tk.Frame(self.root)
        center_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=5, pady=10)
        
        # กราฟ 6 ช่อง (2 แถว 3 คอลัมน์)
        self.fig, self.axs = plt.subplots(2, 3, figsize=(12, 8))
        self.axs = self.axs.flatten()
        self.fig.subplots_adjust(left=0.04, right=0.98, top=0.90, bottom=0.08, hspace=0.45, wspace=0.3)
        
        self.canvas = FigureCanvasTkAgg(self.fig, master=center_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    def draw_scatter(self, ax, title, x_errs, y_errs, tol_m, N, fixed_max_plot):
        theta = np.linspace(0, 2*np.pi, 100)
        ax.plot(tol_m * np.cos(theta), tol_m * np.sin(theta), 'g--', lw=2, label=f'Target ({tol_m*1000:.0f}mm)')
        ax.plot(0, 0, 'k+', markersize=15, lw=2, label='Ideal Center')
        
        hits, miss, max_err = 0, 0, 0
        hit_x, hit_y, miss_x, miss_y = [], [], [], []
        valid_distances = []
        
        for x, y in zip(x_errs, y_errs):
            d = math.hypot(x, y)
            if d >= 5: 
                miss += 1
                continue
            valid_distances.append(d)
            max_err = max(max_err, d)
            if d <= tol_m:
                hits += 1; hit_x.append(x); hit_y.append(y)
            else:
                miss += 1; miss_x.append(x); miss_y.append(y)
                
        hr = (hits/N)*100 if N > 0 else 0
        mean_err = np.mean(valid_distances) if valid_distances else 0
        std_err = np.std(valid_distances) if valid_distances else 0
        
        ax.scatter(hit_x, hit_y, color='green', alpha=0.6, s=15, label='HIT')
        ax.scatter(miss_x, miss_y, color='red', alpha=0.6, s=15, label='MISS')
        
        # ปรับขนาดฟอนต์ Title เพื่อไม่ให้ล้นในช่องแคบ
        ax.set_title(f"{title}\nHit Rate: {hr:.1f}%", fontweight='bold', fontsize=10)
        ax.set_xlabel("Err X (m)", fontsize=9); ax.set_ylabel("Err Y (m)", fontsize=9)
        ax.set_aspect('equal')
        ax.grid(True, linestyle=':', alpha=0.6)
        ax.tick_params(axis='both', which='major', labelsize=8)
        
        ax.set_xlim(-fixed_max_plot, fixed_max_plot)
        ax.set_ylim(-fixed_max_plot, fixed_max_plot)
        ax.legend(loc='upper right', fontsize=7)
        
        return hits, miss, hr, max_err, mean_err, std_err

    def calculate_analytical_bounds(self, v0, p, y, p0, tz, g, R, e_pitch, e_yaw):
        """คำนวณหาสมการและขอบเขตชดเชยระหว่าง Pitch และ Yaw"""
        try:
            res_base = PhysicsEngine.forward_kinematics(v0, p, y, p0, tz, g)
            res_plus = PhysicsEngine.forward_kinematics(v0, p + 0.1, y, p0, tz, g)
            L = res_base['distance']
            S_pitch = abs(res_plus['distance'] - L) / 0.1 # เมตร ต่อ องศา
        except:
            return "Error calculating bounds.", {}

        max_pitch_alone = R / S_pitch if S_pitch > 0 else 0
        max_yaw_alone = math.degrees(math.atan(R / L))

        yaw_dev = L * math.tan(math.radians(e_yaw))
        if yaw_dev >= R:
            allowed_pitch = 0.0
        else:
            allowed_pitch = math.sqrt(R**2 - yaw_dev**2) / S_pitch

        pitch_dev = S_pitch * e_pitch
        if pitch_dev >= R:
            allowed_yaw = 0.0
        else:
            allowed_yaw = math.degrees(math.atan(math.sqrt(R**2 - pitch_dev**2) / L))

        info_text = (
            f"Distance (L) = {L:.2f} m\n"
            f"Target Rad (R) = {R*1000:.0f} mm\n"
            f"Pitch Sensitivity (Sp) = {S_pitch*1000:.1f} mm/deg\n\n"
            f"[Absolute Limits]\n"
            f"Max Pitch (if Yaw=0) = ±{max_pitch_alone:.3f}°\n"
            f"Max Yaw (if Pitch=0) = ±{max_yaw_alone:.3f}°\n\n"
            f"[Conditional Limits]\n"
            f"Given Yaw Err ±{e_yaw}° --> Allowed Pitch ±{allowed_pitch:.3f}°\n"
            f"Given Pitch Err ±{e_pitch}° --> Allowed Yaw ±{allowed_yaw:.3f}°\n\n"
            f"Eq: P_allow = sqrt(R² - (L*tan(Yaw))²) / Sp"
        )
        self.lbl_math.config(text=info_text)
        
        return info_text, {'L': L, 'R': R, 'Sp': S_pitch, 'max_p': max_pitch_alone, 'max_y': max_yaw_alone}

    def run_simulation(self):
        try:
            N = int(self.entries["n_shots"].get())
            tol_m = float(self.entries["tolerance"].get()) / 1000.0
            e_pitch = float(self.entries["err_pitch"].get())
            e_yaw = float(self.entries["err_yaw"].get())
            e_press = float(self.entries["err_press"].get())
            e_delay = float(self.entries["err_delay"].get())
        except ValueError:
            messagebox.showerror("Input Error", "Please enter valid numeric values.")
            return

        for ax in self.axs: ax.clear()
        for item in self.tree_res.get_children(): self.tree_res.delete(item)
        
        base_p, base_y = self.sd['pitch'], self.sd['yaw']
        base_pr, base_v0 = self.sd['pressure'], self.sd['v0']
        p0, g, target_z = self.sd['p0'], self.sd['g'], self.sd['impact_point'][2]
        ix_ideal, iy_ideal = self.sd['impact_point'][0], self.sd['impact_point'][1]
        raw = self.sd['raw_inputs']
        bore, stroke, mass = float(raw['bore']), float(raw['stroke']), float(raw['mass'])
        total_t, math_trig, omega = self.sd['t_total_sync'], self.sd['math_trigger'], self.sd['omega']
        cx, cy, r_tgt = 2.5, 0.0, self.sd['target_radius']

        # คำนวณสมการ Analytical
        math_text, math_data = self.calculate_analytical_bounds(base_v0, base_p, base_y, p0, target_z, g, tol_m, e_pitch, e_yaw)

        # 6 Cases Arrays
        c_x = [[] for _ in range(6)]; c_y = [[] for _ in range(6)]

        for _ in range(N):
            pe = base_p + np.random.uniform(-e_pitch, e_pitch)
            ye = base_y + np.random.uniform(-e_yaw, e_yaw)
            pr_err = max(base_pr + np.random.uniform(-e_press, e_press), 0.01)
            d_err = np.random.uniform(-e_delay, e_delay)

            # Case 1: Pitch Only
            try:
                r1 = PhysicsEngine.forward_kinematics(base_v0, pe, base_y, p0, target_z, g)
                c_x[0].append(r1['impact_point'][0] - ix_ideal); c_y[0].append(r1['impact_point'][1] - iy_ideal)
            except: c_x[0].append(999); c_y[0].append(999)

            # Case 2: Yaw Only
            try:
                r2 = PhysicsEngine.forward_kinematics(base_v0, base_p, ye, p0, target_z, g)
                c_x[1].append(r2['impact_point'][0] - ix_ideal); c_y[1].append(r2['impact_point'][1] - iy_ideal)
            except: c_x[1].append(999); c_y[1].append(999)

            # Case 3: Pitch + Yaw Combined
            try:
                r3 = PhysicsEngine.forward_kinematics(base_v0, pe, ye, p0, target_z, g)
                c_x[2].append(r3['impact_point'][0] - ix_ideal); c_y[2].append(r3['impact_point'][1] - iy_ideal)
            except: c_x[2].append(999); c_y[2].append(999)

            # Case 4: Pressure Only
            try:
                v0_e, _, _ = PhysicsEngine.calculate_v0(pr_err, bore, stroke, mass)
                r4 = PhysicsEngine.forward_kinematics(v0_e, base_p, base_y, p0, target_z, g)
                c_x[3].append(r4['impact_point'][0] - ix_ideal); c_y[3].append(r4['impact_point'][1] - iy_ideal)
            except: c_x[3].append(999); c_y[3].append(999)

            # Case 5: Delay Only
            t_arr = total_t + d_err
            ang5 = math.radians((math_trig + omega * t_arr) % 360)
            c_x[4].append(ix_ideal - (cx + r_tgt * math.cos(ang5)))
            c_y[4].append(iy_ideal - (cy + r_tgt * math.sin(ang5)))

            # Case 6: ALL Combined
            try:
                v0_e, _, _ = PhysicsEngine.calculate_v0(pr_err, bore, stroke, mass)
                r6 = PhysicsEngine.forward_kinematics(v0_e, pe, ye, p0, target_z, g)
                ang6 = math.radians((math_trig + omega * (total_t + d_err))) % 360
                c_x[5].append(r6['impact_point'][0] - (cx + r_tgt * math.cos(ang6)))
                c_y[5].append(r6['impact_point'][1] - (cy + r_tgt * math.sin(ang6)))
            except: c_x[5].append(999); c_y[5].append(999)

        cases_data = [(c_x[i], c_y[i]) for i in range(6)]
        
        # --- กำหนดชื่อกราฟให้แสดงค่า Error ที่เซ็ตไว้ ยกเว้นกราฟสุดท้าย ---
        case_names = [
            f"1: Pitch Err (±{e_pitch}°)", 
            f"2: Yaw Err (±{e_yaw}°)", 
            f"3: Combined Angle (±{e_pitch}°, ±{e_yaw}°)", 
            f"4: Pressure Err (±{e_press} bar)", 
            f"5: Delay Err (±{e_delay} s)", 
            "6: ALL ERRORS COMBINED"
        ]
        file_names = ["Case1_Pitch", "Case2_Yaw", "Case3_CombAngles", "Case4_Pressure", "Case5_Delay", "Case6_ALL"]
        stats = []

        # หาค่า Max เพื่อล็อก Scale
        global_max = 0
        for dx, dy in cases_data:
            for x, y in zip(dx, dy):
                d = math.hypot(x, y)
                if d < 5: global_max = max(global_max, d)
        fixed_max_plot = max(tol_m * 2.5, global_max * 1.1)

        # พล็อต 6 กราฟลงช่อง
        for i in range(6):
            hits, miss, hr, m_err, mean_err, std_err = self.draw_scatter(self.axs[i], case_names[i], cases_data[i][0], cases_data[i][1], tol_m, N, fixed_max_plot)
            stats.append((hits, miss, hr, m_err, mean_err, std_err))
            self.tree_res.insert("", "end", values=(case_names[i], f"{hr:.1f}%", f"{m_err*1000:.1f}", f"{mean_err*1000:.1f}"))

        self.canvas.draw()

        self.last_run_data = {
            'inputs': {'N': N, 'tol': tol_m, 'e_pitch': e_pitch, 'e_yaw': e_yaw, 'e_press': e_press, 'e_delay': e_delay},
            'cases_data': cases_data,
            'case_names': case_names,
            'file_names': file_names,
            'stats': stats,
            'fixed_max_plot': fixed_max_plot,
            'math_text': math_text
        }
        
        self.btn_export.config(state=tk.NORMAL)

    def export_results(self):
        if not self.last_run_data: return
            
        folder_idx = 1
        while True:
            folder_name = f"MC_sim_{folder_idx}"
            if not os.path.exists(folder_name):
                os.makedirs(folder_name)
                break
            folder_idx += 1
            
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        csv_filename = os.path.join(folder_name, f"Statistics_{timestamp}.csv")
        txt_filename = os.path.join(folder_name, f"Summary_Report_{timestamp}.txt")
        
        inp = self.last_run_data['inputs']
        cases_data = self.last_run_data['cases_data']
        case_names = self.last_run_data['case_names']
        file_names = self.last_run_data['file_names']
        stats = self.last_run_data['stats']
        fixed_max_plot = self.last_run_data['fixed_max_plot']
        math_txt = self.last_run_data['math_text']

        try:
            # 1. Export CSV
            with open(csv_filename, 'w', newline='', encoding='utf-8') as f_csv:
                writer = csv.writer(f_csv)
                writer.writerow(["Case", "Total_Shots", "Hits", "Misses", "Hit_Rate_Pct", "Max_Error_mm", "Mean_Error_mm", "Std_Dev_mm"])
                for i in range(6):
                    hits, miss, hr, m_err, mean_err, std_err = stats[i]
                    writer.writerow([case_names[i], inp['N'], hits, miss, f"{hr:.2f}", f"{m_err*1000:.2f}", f"{mean_err*1000:.2f}", f"{std_err*1000:.2f}"])
                    
                    # 2. Export Individual Plots
                    fig_ind, ax_ind = plt.subplots(figsize=(6, 6))
                    self.draw_scatter(ax_ind, case_names[i], cases_data[i][0], cases_data[i][1], inp['tol'], inp['N'], fixed_max_plot)
                    fig_ind.tight_layout()
                    fig_ind.savefig(os.path.join(folder_name, f"Plot_{file_names[i]}_{timestamp}.png"), dpi=300)
                    plt.close(fig_ind)

            # 3. Export Combined Plot (6 กราฟ)
            combined_filename = f"Plot_All_Combined_{timestamp}.png"
            self.fig.savefig(os.path.join(folder_name, combined_filename), dpi=300)

            # 4. Export TXT Summary
            with open(txt_filename, 'w', encoding='utf-8') as f_txt:
                f_txt.write("========================================================\n")
                f_txt.write(" MONTE CARLO SIMULATION - STATISTICAL SUMMARY REPORT\n")
                f_txt.write("========================================================\n")
                f_txt.write(f"Date & Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
                
                f_txt.write("[ BASE PARAMETERS ]\n")
                f_txt.write(f"- Pitch Angle  : {self.sd.get('pitch',0):.2f} deg\n")
                f_txt.write(f"- Yaw Angle    : {self.sd.get('yaw',0):.2f} deg\n")
                f_txt.write(f"- Air Pressure : {self.sd.get('pressure',0)} bar\n")
                f_txt.write(f"- Firing Delay : {self.sd.get('t_delay',0):.3f} sec\n\n")

                f_txt.write("[ ANALYTICAL BOUNDARIES ]\n")
                f_txt.write(f"{math_txt}\n\n")

                f_txt.write("[ ERROR SETTINGS (NOISE) ]\n")
                f_txt.write(f"- Total Shots    : {inp['N']} shots per case\n")
                f_txt.write(f"- Target Radius  : {inp['tol']*1000:.0f} mm\n")
                f_txt.write(f"- Pitch Error    : +/- {inp['e_pitch']} deg\n")
                f_txt.write(f"- Yaw Error      : +/- {inp['e_yaw']} deg\n")
                f_txt.write(f"- Pressure Error : +/- {inp['e_press']} bar\n")
                f_txt.write(f"- Delay Error    : +/- {inp['e_delay']} sec\n\n")

                f_txt.write("[ STATISTICAL RESULTS ]\n")
                for i in range(6):
                    hits, miss, hr, m_err, mean_err, std_err = stats[i]
                    f_txt.write(f"--- {case_names[i]} ---\n")
                    f_txt.write(f"  > Hits / Misses : {hits} / {miss}\n")
                    f_txt.write(f"  > Hit Rate      : {hr:.2f} %\n")
                    f_txt.write(f"  > Max Error     : {m_err*1000:.2f} mm\n")
                    f_txt.write(f"  > Mean Error    : {mean_err*1000:.2f} mm\n")
                    f_txt.write(f"  > Std Deviation : {std_err*1000:.2f} mm\n\n")

                f_txt.write("================ END OF REPORT ================\n")

            # 5. Success Message
            msg = (f"Export completed successfully!\n\n"
                   f"Saved 9 files in folder: '{folder_name}'\n"
                   f"📄 Statistics_{timestamp}.csv\n"
                   f"📄 Summary_Report_{timestamp}.txt\n"
                   f"🖼️ 6 Individual Plots (.png)\n"
                   f"🖼️ {combined_filename}")
            messagebox.showinfo("Export Success", msg)
            
        except Exception as e:
            messagebox.showerror("Export Error", f"An error occurred during export:\n{e}")

if __name__ == "__main__":
    app = RepeatabilityWindow()
    app.root.mainloop()