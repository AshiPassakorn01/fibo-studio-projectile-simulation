import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import subprocess
import sys
import json
import tools
from physic import PhysicsEngine, TimingLogic

class ControlWindow:
    def __init__(self):
        self.root = tk.Tk()
        self.root.title("Projectile & Timing Simulator")
        self.root.geometry("480x800") 
        
        self.data = tools.load_data()
        self.mode_var = tk.StringVar(value=self.data.get('mode', 'inverse'))
        self.entries = {}
        self.setup_ui()
        
    def setup_ui(self):
        mode_frame = tk.LabelFrame(self.root, text=" 1. Mode Selection ", font=("Arial", 10, "bold"), padx=10, pady=5)
        mode_frame.pack(fill=tk.X, padx=15, pady=10)
        
        tk.Radiobutton(mode_frame, text="Inverse Mode (Calculate Angle from Target)", variable=self.mode_var, value="inverse", font=("Arial", 9), command=self.update_table).grid(row=0, column=0, sticky='w', pady=2)
        tk.Radiobutton(mode_frame, text="Forward Mode (Calculate Target from Angle)", variable=self.mode_var, value="forward", font=("Arial", 9), command=self.update_table).grid(row=1, column=0, sticky='w', pady=2)

        self.param_frame = tk.LabelFrame(self.root, text=" 2. System Parameters ", font=("Arial", 10, "bold"), padx=10, pady=5)
        self.param_frame.pack(fill=tk.BOTH, expand=True, padx=15, pady=0)
        
        self.canvas = tk.Canvas(self.param_frame, highlightthickness=0)
        scrollbar = ttk.Scrollbar(self.param_frame, orient="vertical", command=self.canvas.yview)
        self.table_frame = tk.Frame(self.canvas)
        
        self.table_frame.bind("<Configure>", lambda e: self.canvas.configure(scrollregion=self.canvas.bbox("all")))
        self.canvas.create_window((0, 0), window=self.table_frame, anchor="nw")
        self.canvas.configure(yscrollcommand=scrollbar.set)
        
        self.canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")

        btn_frame = tk.Frame(self.root)
        btn_frame.pack(fill=tk.X, pady=10)
        
        tk.Button(btn_frame, text="🚀 Run Simulation", font=("Arial", 12, "bold"), bg="#4CAF50", fg="white", cursor="hand2", command=self.run_simulation, width=22, height=2).pack(pady=5)
        
        file_btn_frame = tk.Frame(btn_frame)
        file_btn_frame.pack(pady=5)
        tk.Button(file_btn_frame, text="💾 Save Result", font=("Arial", 10, "bold"), bg="#0275d8", fg="white", cursor="hand2", command=self.save_simulation, width=12).grid(row=0, column=0, padx=5)
        tk.Button(file_btn_frame, text="📂 Load Sim", font=("Arial", 10, "bold"), bg="#f0ad4e", fg="white", cursor="hand2", command=self.load_simulation, width=12).grid(row=0, column=1, padx=5)

        self.update_table()

    def update_table(self):
        for widget in self.table_frame.winfo_children(): widget.destroy()
        self.entries.clear()
        
        tk.Label(self.table_frame, text="Parameter Name", font=("Arial", 9, "bold")).grid(row=0, column=0, pady=(0, 5), sticky='w', padx=5)
        tk.Label(self.table_frame, text="Value", font=("Arial", 9, "bold")).grid(row=0, column=1, pady=(0, 5), sticky='w', padx=5)
        tk.Label(self.table_frame, text="Unit", font=("Arial", 9, "bold")).grid(row=0, column=2, pady=(0, 5), sticky='w', padx=5)

        common_params = [
            ("g", "Gravity Accel. (g):", "9.80665", "m/s²"),
            ("pressure", "Air Pressure:", "0.55", "bar"),
            ("bore", "Cylinder Bore:", "25", "mm"),
            ("stroke", "Cylinder Stroke:", "100", "mm"),
            ("mass", "Total Mass:", "0.150", "kg"),
            ("separator1", "", "", ""), 
            ("x0", "Muzzle Pos X0:", "0.0", "m"),
            ("y0", "Muzzle Pos Y0:", "0.0", "m"),
            ("z0", "Muzzle Pos Z0:", "0.0", "m"),
            ("separator2", "", "", ""),
            ("target_dir", "Rotation Dir (CCW/CW):", "CCW", "-"),
            ("target_rpm", "Target Speed (RPM):", "30", "rpm"),
            ("target_radius", "Target Radius:", "0.15", "m"),
            ("trigger_angle", "Trigger (from closest, CCW):", "0", "deg"),
            ("hw_latency", "Hardware Latency:", "0.15", "s"),
            ("tolerance", "Hit Tolerance:", "5", "deg"),
            ("separator3", "", "", "")  
        ]
        
        inverse_params = [("tz", "Target Z-Level:", "0.0", "m")]
        forward_params = [("pitch", "Pitch Angle:", "45.0", "deg"), ("yaw", "Yaw Angle:", "0.0", "deg"), ("target_z", "Impact Z-Level:", "0.0", "m")]

        current_mode = self.mode_var.get()
        active_params = common_params.copy() + (inverse_params if current_mode == "inverse" else forward_params)

        for i, (key, name, default, unit) in enumerate(active_params, start=1):
            if key.startswith("separator"):
                ttk.Separator(self.table_frame, orient='horizontal').grid(row=i, column=0, columnspan=3, sticky='ew', pady=5)
                continue
            tk.Label(self.table_frame, text=name, font=("Arial", 9)).grid(row=i, column=0, sticky='w', padx=5, pady=2)
            ent = tk.Entry(self.table_frame, width=10, font=("Consolas", 10), justify='left')
            ent.insert(0, default)
            ent.grid(row=i, column=1, padx=5, pady=2, sticky='w')
            tk.Label(self.table_frame, text=unit, font=("Arial", 9), fg="#555555").grid(row=i, column=2, sticky='w', padx=5, pady=2)
            self.entries[key] = ent

        if 'raw_inputs' in self.data:
            for k, v in self.data['raw_inputs'].items():
                if k in self.entries:
                    self.entries[k].delete(0, tk.END)
                    self.entries[k].insert(0, str(v))

    def _calculate_sim_data(self):
        d = {k: v.get() if k == "target_dir" else float(v.get()) for k, v in self.entries.items()}
        current_mode = self.mode_var.get()
        
        v0, work, force = PhysicsEngine.calculate_v0(d['pressure'], d['bore'], d['stroke'], d['mass'])
        p0, g = [d['x0'], d['y0'], d['z0']], d['g']
        
        math_impact = 180.0
        math_trigger = (180.0 + d['trigger_angle']) % 360.0
        
        sim_data = {
            'mode': current_mode, 'pressure': d['pressure'],
            'v0': v0, 'work': work, 'force': force, 'p0': p0, 'g': g,
            'target_dir': str(d['target_dir']).strip().upper(),
            'target_rpm': d['target_rpm'], 'target_radius': d['target_radius'],
            'trigger_angle': d['trigger_angle'], 'math_trigger': math_trigger, 'math_impact': math_impact,
            'hw_latency': d['hw_latency'], 'tolerance': d['tolerance'], 'raw_inputs': d 
        }
        
        if current_mode == 'inverse':
            p_target = [2.5 - d['target_radius'], 0.0, d['tz']]
            res = PhysicsEngine.inverse_kinematics(v0, p0, p_target, g)
            sim_data.update({'pitch': res['pitch'], 'yaw': res['yaw'], 't_flight': res['t_flight'], 'distance': res['distance'], 'impact_point': p_target})
        else:
            res = PhysicsEngine.forward_kinematics(v0, d['pitch'], d['yaw'], p0, d['target_z'], g)
            sim_data.update({'pitch': d['pitch'], 'yaw': d['yaw'], 't_flight': res['t_flight'], 'distance': res['distance'], 'impact_point': res['impact_point']})
            
        t_delay, t_total, omega = TimingLogic.calculate_delay(d['target_rpm'], math_trigger, math_impact, d['hw_latency'], sim_data['t_flight'], sim_data['target_dir'])
        sim_data.update({'t_delay': t_delay, 't_total_sync': t_total, 'omega': omega})
            
        return sim_data

    def run_simulation(self):
        try:
            tools.update_data(self._calculate_sim_data())
            subprocess.Popen([sys.executable, "win3.py"])
        except Exception as e:
            messagebox.showerror("Error", f"Error:\n{e}")

    def save_simulation(self):
        try:
            filepath = filedialog.asksaveasfilename(defaultextension=".json", filetypes=[("JSON Files", "*.json")])
            if filepath:
                with open(filepath, 'w', encoding='utf-8') as f: json.dump(self._calculate_sim_data(), f, indent=4)
                messagebox.showinfo("Success", "Saved successfully!")
        except Exception as e: messagebox.showerror("Error", str(e))

    def load_simulation(self):
        filepath = filedialog.askopenfilename(filetypes=[("JSON Files", "*.json")])
        if filepath:
            try:
                with open(filepath, 'r', encoding='utf-8') as f: loaded_data = json.load(f)
                tools.update_data(loaded_data)
                self.data = loaded_data
                if 'mode' in loaded_data:
                    self.mode_var.set(loaded_data['mode'])
                    self.update_table()
                subprocess.Popen([sys.executable, "win3.py"])
            except Exception as e: messagebox.showerror("Load Error", str(e))

    def run(self): self.root.mainloop()

if __name__ == "__main__":
    app = ControlWindow()
    app.run()