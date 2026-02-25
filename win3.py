import tkinter as tk
from tkinter import ttk
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.animation import FuncAnimation
import math
import tools
from physic import TimingLogic

class SimulationWindow:
    def __init__(self):
        self.root = tk.Tk()
        self.root.title("3D Trajectory & Timing Sync Animator (Smooth Request Mode)")
        self.root.geometry("1450x980") 
        self.sd = tools.load_data()
        
        # --- ลอจิกควบคุมการยิงและสัญญาณ ---
        self.is_playing = False
        self.request_armed = False        # สถานะการกดปุ่มขอสัญญาณค้างไว้
        self.firing_active = False       # สถานะการเริ่มนับเวลาและยิงลูกบอล
        self.fire_frame_counter = 0      # ตัวนับเฟรมของลูกบอล (เพื่อดึงค่าล่วงหน้าที่คำนวณไว้)
        self.last_trigger_time = -1.0    # ป้องกันการส่งสัญญาณซ้ำซ้อนในรอบเดียวกัน
        
        self.setup_ui()
        self.setup_plots()
        
    def setup_ui(self):
        info_frame = tk.Frame(self.root, width=460, bg="#f0f0f0", relief=tk.RIDGE, bd=2)
        info_frame.pack(side=tk.LEFT, fill=tk.Y, padx=10, pady=10)
        info_frame.pack_propagate(False) 
        
        # 1. Inputs & 2. Results
        for i, title in enumerate([" 1. Initial Parameters ", " 2. Calculated Results "]):
            frame = tk.LabelFrame(info_frame, text=title, font=("Arial", 10, "bold"), bg="#f0f0f0")
            frame.pack(fill=tk.X, padx=10, pady=2)
            tree = ttk.Treeview(frame, columns=("Param", "Value"), show="headings", height=5)
            tree.heading("Param", text="Parameter"); tree.column("Param", width=180)
            tree.heading("Value", text="Value"); tree.column("Value", width=200)
            tree.pack(fill=tk.X, padx=5, pady=2)
            if i == 0: self.tree_in = tree
            else: self.tree_out = tree

        # 3. Position Log (5 บรรทัด)
        log_frame = tk.LabelFrame(info_frame, text=" 3. Position Log (Δt = 0.05s) ", font=("Arial", 10, "bold"), bg="#f0f0f0")
        log_frame.pack(fill=tk.X, padx=10, pady=2)
        self.tree_log = ttk.Treeview(log_frame, columns=("Time", "X", "Y", "Z"), show="headings", height=5)
        for col in ("Time", "X", "Y", "Z"):
            self.tree_log.heading(col, text=col); self.tree_log.column(col, width=75, anchor="center")
        self.tree_log.pack(fill=tk.X, padx=5, pady=2)

        # 4. Real-time Kinematics Monitor
        mon_frame = tk.LabelFrame(info_frame, text=" 4. Real-time Kinematics Monitor ", font=("Arial", 10, "bold"), bg="#e8f4f8", fg="#005b96")
        mon_frame.pack(fill=tk.X, padx=10, pady=2)
        self.lbl_mon_time = tk.Label(mon_frame, text="Time: 0.000 s | State: -", font=("Consolas", 10, "bold"), bg="#e8f4f8")
        self.lbl_mon_time.pack(anchor="w", padx=10, pady=2)
        self.tree_mon = ttk.Treeview(mon_frame, columns=("Axis", "Pos", "Vel", "Acc"), show="headings", height=3)
        for col in ("Axis", "Pos", "Vel", "Acc"):
            self.tree_mon.heading(col, text=col); self.tree_mon.column(col, width=90, anchor="center")
        self.tree_mon.pack(fill=tk.X, padx=5, pady=2)
        for ax in ["X", "Y", "Z"]: self.tree_mon.insert("", "end", iid=ax, values=(ax, "0.000", "0.000", "0.000"))

        self.populate_tables()

        # 5. View & System Controls
        view_frame = tk.LabelFrame(info_frame, text=" 5. System Controls ", bg="#f0f0f0", font=("Arial", 10, "bold"))
        view_frame.pack(fill=tk.X, padx=10, pady=5)
        
        # --- ส่วนของปุ่ม Request (ต้องกดค้างเพื่อขอสัญญาณ) ---
        self.btn_req = tk.Button(view_frame, text="🔘 HOLD TO REQUEST TRIGGER", font=("Arial", 10, "bold"), bg="#f0ad4e", fg="white", height=2)
        self.btn_req.pack(fill=tk.X, padx=20, pady=5)
        self.btn_req.bind("<ButtonPress-1>", self.on_req_press)
        self.btn_req.bind("<ButtonRelease-1>", self.on_req_release)

        v_btn_frame = tk.Frame(view_frame, bg="#f0f0f0")
        v_btn_frame.pack(pady=2)
        tk.Button(v_btn_frame, text="Top View", width=12,fg="black", bg="#a0ffa9" , command=lambda: self.change_view(90, -90)).grid(row=0, column=0, padx=2)
        tk.Button(v_btn_frame, text="Front View", width=12,fg="black", bg="#ff7b7b" , command=lambda: self.change_view(0, -90)).grid(row=0, column=1, padx=2)
        tk.Button(v_btn_frame, text="Side View", width=12,fg="black", bg="#a0b4ff", command=lambda: self.change_view(0, 0)).grid(row=0, column=2, padx=2)
        tk.Button(view_frame, text="Default 3D View", width=40, command=lambda: self.change_view(30, -60)).pack(pady=2)

        ctrl_frame = tk.Frame(view_frame, bg="#f0f0f0")
        ctrl_frame.pack(fill=tk.X, pady=5)
        tk.Button(ctrl_frame, text="Start Simulation", font=("Arial", 10, "bold"), bg="#5cb85c", fg="white", command=self.play_anim).pack(side=tk.TOP, fill=tk.X, padx=20, pady=2)
        tk.Button(ctrl_frame, text="Stop", font=("Arial", 10, "bold"), bg="#888", fg="white", command=self.pause_anim).pack(side=tk.TOP, fill=tk.X, padx=20, pady=2)
        tk.Button(ctrl_frame, text="Reset System", font=("Arial", 10, "bold"), bg="#d9534f", fg="white", command=self.reset_anim).pack(side=tk.TOP, fill=tk.X, padx=20, pady=2)

        self.right_frame = tk.Frame(self.root)
        self.right_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True, padx=10, pady=10)
        self.notebook = ttk.Notebook(self.right_frame)
        self.notebook.pack(fill=tk.BOTH, expand=True)
        self.tab_3d = tk.Frame(self.notebook); self.tab_sync = tk.Frame(self.notebook)
        self.notebook.add(self.tab_3d, text="3D Trajectory "); self.notebook.add(self.tab_sync, text=" 🎯 3D Timing Sync ")

    def on_req_press(self, e): 
        self.request_armed = True
        self.btn_req.config(bg="#28a745", text="REQUEST ACTIVE (WAITING TRIGGER)")

    def on_req_release(self, e): 
        self.request_armed = False
        if not self.firing_active:
            self.btn_req.config(bg="#f0ad4e", text="REQUEST TRIGGER")

    def populate_tables(self):
        inputs = [("Mode", self.sd.get('mode','').upper()), ("Pressure", f"{self.sd.get('pressure')} bar"), ("Gravity", f"{self.sd.get('g'):.4f}"), ("RPM", f"{self.sd.get('target_rpm')} ({self.sd.get('target_dir')})"), ("HW Latency", f"{self.sd.get('hw_latency')} s")]
        calcs = [("V0", f"{self.sd.get('v0',0):.3f} m/s"), ("T-Flight", f"{self.sd.get('t_flight',0):.3f} s"), ("T-Delay", f"{self.sd.get('t_delay',0):.3f} s"), ("Pitch", f"{self.sd.get('pitch',0):.2f}°"), ("Yaw", f"{self.sd.get('yaw',0):.2f}°")]
        for item in inputs: self.tree_in.insert("", "end", values=item)
        for item in calcs: self.tree_out.insert("", "end", values=item)
        t_f = self.sd['t_flight']
        t_steps = np.arange(0, t_f, 0.05)
        if len(t_steps)==0 or t_steps[-1]<t_f: t_steps = np.append(t_steps, t_f)
        v0, pr, yr, p0 = self.sd['v0'], np.radians(self.sd['pitch']), np.radians(self.sd['yaw']), self.sd['p0']
        for t in t_steps:
            x = p0[0] + (v0*np.cos(pr)*np.cos(yr))*t
            y = p0[1] + (v0*np.cos(pr)*np.sin(yr))*t
            z = max(p0[2] + (v0*np.sin(pr))*t - 0.5*self.sd['g']*t**2, self.sd['impact_point'][2] if t>=t_f else -999)
            self.tree_log.insert("", "end", values=(f"{t:.3f}", f"{x:.3f}", f"{y:.3f}", f"{z:.3f}"))

    def style_3d_axes(self, ax):
        ax.set_xlabel("X (m)", color='blue'); ax.set_ylabel("Y (m)", color='red'); ax.set_zlabel("Z (m)", color='green')
        ax.zaxis.pane.set_facecolor((0,1,0,0.05)); ax.yaxis.pane.set_facecolor((1,0,0,0.05)); ax.xaxis.pane.set_facecolor((0,0,1,0.05))
        ax.set_autoscale_on(False)

    def set_axes_limits(self, ax, all_x, all_y, all_z):
        x_m, x_M = min(all_x), max(all_x); y_m, y_M = min(all_y), max(all_y); z_m, z_M = min(all_z), max(all_z)
        rng = max(x_M-x_m, y_M-y_m, z_M-z_m, 1); half = (rng + rng*0.6)*0.5
        mx, my, mz = (x_M+x_m)/2, (y_M+y_m)/2, (z_M+z_m)/2
        ax.set_xlim(mx-half, mx+half); ax.set_ylim(my-half, my+half); ax.set_zlim(0, mz+half)
        ax.set_box_aspect((1,1,1))

    def setup_plots(self):
        # --- Pre-calculate Trajectory Data ---
        t_f = self.sd['t_flight']
        self.t_arr = np.arange(0, t_f, 0.05)
        if len(self.t_arr)==0 or self.t_arr[-1]<t_f: self.t_arr = np.append(self.t_arr, t_f)
        v0, g, pr, yr, p0 = self.sd['v0'], self.sd['g'], np.radians(self.sd['pitch']), np.radians(self.sd['yaw']), self.sd['p0']
        self.x_3d = p0[0] + (v0*np.cos(pr)*np.cos(yr))*self.t_arr
        self.y_3d = p0[1] + (v0*np.cos(pr)*np.sin(yr))*self.t_arr
        self.z_3d = p0[2] + (v0*np.sin(pr))*self.t_arr - 0.5*g*self.t_arr**2

        # --- Pre-calculate Sync Animation Data (Smooth Projectile) ---
        self.fps = 20
        self.total_t_sync = self.sd['t_total_sync']
        self.frames_projectile = int(np.ceil(self.total_t_sync * self.fps)) + 1
        
        self.sync_bx = np.zeros(self.frames_projectile); self.sync_by = np.zeros(self.frames_projectile); self.sync_bz = np.zeros(self.frames_projectile)
        self.sync_b_vis = np.zeros(self.frames_projectile, dtype=bool)
        self.sync_titles, self.sync_colors = [], []
        t_d, t_h, impact = self.sd['t_delay'], self.sd['hw_latency'], self.sd['impact_point']
        
        for f in range(self.frames_projectile):
            t = f/self.fps
            if t < t_d:
                self.sync_b_vis[f] = False; self.sync_titles.append(f"DELAYING: {t:.2f}s"); self.sync_colors.append("orange")
            elif t < t_d+t_h:
                self.sync_b_vis[f] = True; self.sync_bx[f],self.sync_by[f],self.sync_bz[f] = p0; self.sync_titles.append(f"FIRING (Valve Open): {t:.2f}s"); self.sync_colors.append("red")
            elif t <= self.total_t_sync:
                self.sync_b_vis[f] = True; tp = min(t-(t_d+t_h), t_f)
                self.sync_bx[f] = p0[0] + (v0*np.cos(pr)*np.cos(yr))*tp
                self.sync_by[f] = p0[1] + (v0*np.cos(pr)*np.sin(yr))*tp
                self.sync_bz[f] = max(p0[2] + (v0*np.sin(pr))*tp - 0.5*g*tp**2, impact[2])
                self.sync_titles.append(f"FLYING: {t:.2f}s"); self.sync_colors.append("blue")
            else:
                self.sync_b_vis[f] = True; self.sync_bx[f],self.sync_by[f],self.sync_bz[f] = impact
                msg, color, err = TimingLogic.validate_impact(self.sd['math_trigger'] + self.sd['omega']*self.total_t_sync, self.sd['math_impact'], self.sd['tolerance'])
                self.sync_titles.append(f"IMPACT! [{msg}]"); self.sync_colors.append(color)

        # --- Setup Matplotlib Figures ---
        self.fig_3d = plt.figure(); self.ax_3d = self.fig_3d.add_subplot(111, projection='3d')
        self.canvas_3d = FigureCanvasTkAgg(self.fig_3d, self.tab_3d); self.canvas_3d.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        self.style_3d_axes(self.ax_3d)
        self.line_3d, = self.ax_3d.plot([],[],[], 'b-', lw=2); self.point_3d, = self.ax_3d.plot([],[],[], 'ro')
        self.ax_3d.scatter(*p0, color='g', s=50); self.ax_3d.scatter(*impact, color='r', marker='x', s=100)

        self.fig_sync = plt.figure(); self.ax_sync = self.fig_sync.add_subplot(111, projection='3d')
        self.canvas_sync = FigureCanvasTkAgg(self.fig_sync, self.tab_sync); self.canvas_sync.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        self.style_3d_axes(self.ax_sync)
        self.status_text = self.ax_sync.text2D(0.02, 0.95, "IDLE", transform=self.ax_sync.transAxes, weight='bold', bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))
        
        cx, cy, cz, r = 2.5, 0, impact[2], self.sd['target_radius']
        th = np.linspace(0, 2*np.pi, 100)
        self.ax_sync.plot(cx+r*np.cos(th), cy+r*np.sin(th), np.full_like(th, cz), 'gray', ls='--')
        self.ax_sync.plot([cx], [cy], [cz], 'k^', ms=10)
        trig_rad = np.radians(self.sd['math_trigger'])
        self.ax_sync.plot([cx, cx+r*np.cos(trig_rad)], [cy, cy+r*np.sin(trig_rad)], [cz, cz], 'orange', lw=2)
        
        self.target_sync, = self.ax_sync.plot([],[],[], 'bo', ms=10); self.ball_sync, = self.ax_sync.plot([],[],[], 'go', ms=8)
        
        all_x = np.concatenate([self.x_3d, [cx+r, cx-r]]); all_y = np.concatenate([self.y_3d, [cy+r, cy-r]]); all_z = np.concatenate([self.z_3d, [cz]])
        self.set_axes_limits(self.ax_3d, all_x, all_y, all_z); self.set_axes_limits(self.ax_sync, all_x, all_y, all_z)
        self.change_view(30, -60)
        
        self.anim_3d = FuncAnimation(self.fig_3d, self.update_3d, frames=len(self.t_arr), interval=50)
        self.anim_sync = FuncAnimation(self.fig_sync, self.update_sync, interval=50, cache_frame_data=False)
        self.pause_anim()

    def update_monitor(self, t, tp, state):
        self.lbl_mon_time.config(text=f"Time: {t:.3f} s | State: {state}")
        v0, g, pr, yr, p0 = self.sd['v0'], self.sd['g'], np.radians(self.sd['pitch']), np.radians(self.sd['yaw']), self.sd['p0']
        vx0, vy0, vz0 = v0*np.cos(pr)*np.cos(yr), v0*np.cos(pr)*np.sin(yr), v0*np.sin(pr)
        if "FLYING" in state:
            pos = [p0[0]+vx0*tp, p0[1]+vy0*tp, p0[2]+vz0*tp-0.5*g*tp**2]
            vel = [vx0, vy0, vz0-g*tp]; acc = [0, 0, -g]
        elif "IMPACT" in state:
            pos = self.sd['impact_point']; vel = [0,0,0]; acc = [0,0,0]
        else: pos = p0; vel = [0,0,0]; acc = [0,0,0]
        for i, ax in enumerate(["X", "Y", "Z"]): self.tree_mon.item(ax, values=(ax, f"{pos[i]:.3f}", f"{vel[i]:.3f}", f"{acc[i]:.3f}"))

    def update_3d(self, i):
        self.line_3d.set_data(self.x_3d[:i+1], self.y_3d[:i+1]); self.line_3d.set_3d_properties(self.z_3d[:i+1])
        self.point_3d.set_data([self.x_3d[i]], [self.y_3d[i]]); self.point_3d.set_3d_properties([self.z_3d[i]])
        if self.notebook.index(self.notebook.select())==0: 
            self.update_monitor(self.t_arr[i], self.t_arr[i], "FLYING" if self.t_arr[i]<self.sd['t_flight'] else "IMPACT")
        return self.line_3d, self.point_3d

    def update_sync(self, frame):
        if not self.is_playing: return self.target_sync, self.ball_sync
        
        t_global = frame * 0.05
        cx, cy, cz, r = 2.5, 0, self.sd['impact_point'][2], self.sd['target_radius']
        
        # 1. เป้าหมายหมุนอย่างต่อเนื่อง
        curr_ang = (self.sd['math_trigger'] + (self.sd['omega'] * t_global)) % 360
        self.target_sync.set_data([cx+r*math.cos(math.radians(curr_ang))], [cy+r*math.sin(math.radians(curr_ang))])
        self.target_sync.set_3d_properties([cz])

        # 2. ตรวจจับสัญญาณ Trigger อย่างแม่นยำ (คำนวณจากเวลาครบรอบ T_rev)
        T_rev = 360.0 / abs(self.sd['omega'])
        time_since_trigger_pass = t_global % T_rev
        
        # ถ้าเป้าเพิ่งผ่านจุดทริกเกอร์มาไม่เกิน 1 เฟรม (0.06 วิ) และกดปุ่มค้างอยู่ ให้ปล่อยสัญญาณยิง!
        if time_since_trigger_pass <= 0.06 and self.request_armed and not self.firing_active:
            # ป้องกันการส่งสัญญาณซ้ำรัวๆ ในรอบเดียวกัน
            if t_global - self.last_trigger_time > (T_rev / 2):
                self.firing_active = True
                self.fire_frame_counter = 0
                self.last_trigger_time = t_global
                self.btn_req.config(text="SIGNAL SENT! (FIRING)")

        # ทำเอฟเฟกต์กระพริบสีทองเมื่อเป้าทับเส้น Trigger 
        if time_since_trigger_pass <= 0.15:
            self.target_sync.set_color('#FFD700'); self.target_sync.set_markersize(15)
        else:
            self.target_sync.set_color('blue'); self.target_sync.set_markersize(10)

        # 3. จัดการแอนิเมชันของลูกบอล (ดึงข้อมูลที่ Pre-calculate ไว้มาเล่น เพื่อความลื่นไหล)
        if self.firing_active:
            f = self.fire_frame_counter
            if f < self.frames_projectile:
                if self.sync_b_vis[f]:
                    self.ball_sync.set_data([self.sync_bx[f]], [self.sync_by[f]])
                    self.ball_sync.set_3d_properties([self.sync_bz[f]])
                else:
                    self.ball_sync.set_data([],[]); self.ball_sync.set_3d_properties([])
                
                self.status_text.set_text(self.sync_titles[f]); self.status_text.set_color(self.sync_colors[f])
                
                # อัปเดตตารางมอนิเตอร์เฉพาะตอนที่ดูหน้า Sync
                if self.notebook.index(self.notebook.select())==1:
                    t_p = max(0, f/self.fps - (self.sd['t_delay']+self.sd['hw_latency']))
                    self.update_monitor(f/self.fps, t_p, self.sync_titles[f].split(":")[0])
                
                self.fire_frame_counter += 1
            else:
                # ยิงจบแล้ว รีเซ็ตสถานะกลับไปรอ Request ใหม่
                self.firing_active = False 
                if self.request_armed: self.btn_req.config(text="REQUEST ACTIVE (WAITING TRIGGER)")
        else:
            self.ball_sync.set_data([],[]); self.ball_sync.set_3d_properties([])
            if self.request_armed:
                self.status_text.set_text("ARMED: WAITING FOR TARGET..."); self.status_text.set_color("green")
            else:
                self.status_text.set_text("IDLE: HOLD REQUEST TO ARM"); self.status_text.set_color("black")
        
        return self.target_sync, self.ball_sync

    def change_view(self, e, a):
        self.ax_3d.view_init(elev=e, azim=a); self.ax_sync.view_init(elev=e, azim=a)
        self.canvas_3d.draw(); self.canvas_sync.draw()

    def play_anim(self): self.is_playing = True
    def pause_anim(self): self.is_playing = False
    def reset_anim(self):
        self.is_playing = False; self.firing_active = False; self.fire_frame_counter = 0; self.last_trigger_time = -1.0
        self.btn_req.config(bg="#f0ad4e", text=" REQUEST TRIGGER")
        self.status_text.set_text("SYSTEM RESET"); self.canvas_sync.draw()

if __name__ == "__main__":
    SimulationWindow().root.mainloop()