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
        self.root.title("3D Trajectory & Timing Sync Animator (Perfect Accuracy)")
        self.root.geometry("1450x980") 
        self.sd = tools.load_data()
        self.is_playing = False
        self.setup_ui()
        self.setup_plots()
        
    def setup_ui(self):
        info_frame = tk.Frame(self.root, width=460, bg="#f0f0f0", relief=tk.RIDGE, bd=2)
        info_frame.pack(side=tk.LEFT, fill=tk.Y, padx=10, pady=10)
        info_frame.pack_propagate(False) 
        
        # 1. Inputs
        input_frame = tk.LabelFrame(info_frame, text=" 1. Initial Parameters ", font=("Arial", 10, "bold"), bg="#f0f0f0")
        input_frame.pack(fill=tk.X, padx=10, pady=2)
        self.tree_in = ttk.Treeview(input_frame, columns=("Param", "Value"), show="headings", height=5)
        self.tree_in.heading("Param", text="Parameter"); self.tree_in.column("Param", width=180)
        self.tree_in.heading("Value", text="Value"); self.tree_in.column("Value", width=200)
        self.tree_in.pack(fill=tk.X, padx=5, pady=2)

        # 2. Results
        calc_frame = tk.LabelFrame(info_frame, text=" 2. Calculated Results ", font=("Arial", 10, "bold"), bg="#f0f0f0")
        calc_frame.pack(fill=tk.X, padx=10, pady=2)
        self.tree_out = ttk.Treeview(calc_frame, columns=("Param", "Value"), show="headings", height=5)
        self.tree_out.heading("Param", text="Parameter"); self.tree_out.column("Param", width=180)
        self.tree_out.heading("Value", text="Value"); self.tree_out.column("Value", width=200)
        self.tree_out.pack(fill=tk.X, padx=5, pady=2)

        # 3. Position Log 
        log_frame = tk.LabelFrame(info_frame, text=" 3. Position Log (Δt = 0.05s) ", font=("Arial", 10, "bold"), bg="#f0f0f0")
        log_frame.pack(fill=tk.X, padx=10, pady=2)
        self.tree_log = ttk.Treeview(log_frame, columns=("Time", "X", "Y", "Z"), show="headings", height=5)
        for col in ("Time", "X", "Y", "Z"):
            self.tree_log.heading(col, text=f"{col} {'(s)' if col=='Time' else '(m)'}")
            self.tree_log.column(col, width=75, anchor="center")
        scroll = ttk.Scrollbar(log_frame, orient=tk.VERTICAL, command=self.tree_log.yview)
        self.tree_log.configure(yscroll=scroll.set)
        scroll.pack(side=tk.RIGHT, fill=tk.Y)
        self.tree_log.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, pady=2)
        self.populate_tables()

        # 4. Real-time Kinematics Monitor 
        mon_frame = tk.LabelFrame(info_frame, text=" 4. Real-time Kinematics Monitor ", font=("Arial", 10, "bold"), bg="#e8f4f8", fg="#005b96")
        mon_frame.pack(fill=tk.X, padx=10, pady=2)
        
        self.lbl_mon_time = tk.Label(mon_frame, text="Time: 0.000 s | State: -", font=("Consolas", 10, "bold"), bg="#e8f4f8", fg="#000")
        self.lbl_mon_time.pack(anchor="w", padx=10, pady=2)
        
        self.tree_mon = ttk.Treeview(mon_frame, columns=("Axis", "Pos", "Vel", "Acc"), show="headings", height=3)
        self.tree_mon.heading("Axis", text="Axis")
        self.tree_mon.heading("Pos", text="Pos (m)")
        self.tree_mon.heading("Vel", text="Vel (m/s)")
        self.tree_mon.heading("Acc", text="Acc (m/s²)")
        self.tree_mon.column("Axis", width=40, anchor="center")
        self.tree_mon.column("Pos", width=90, anchor="center")
        self.tree_mon.column("Vel", width=90, anchor="center")
        self.tree_mon.column("Acc", width=90, anchor="center")
        self.tree_mon.pack(fill=tk.X, padx=5, pady=2)
        
        self.tree_mon.insert("", "end", iid="X", values=("X", "0.000", "0.000", "0.000"))
        self.tree_mon.insert("", "end", iid="Y", values=("Y", "0.000", "0.000", "0.000"))
        self.tree_mon.insert("", "end", iid="Z", values=("Z", "0.000", "0.000", "0.000"))

        # 5. View & Controls
        view_frame = tk.LabelFrame(info_frame, text=" 5. View & System Controls ", bg="#f0f0f0", font=("Arial", 10, "bold"))
        view_frame.pack(fill=tk.X, padx=10, pady=5)
        v_btn_frame1 = tk.Frame(view_frame, bg="#f0f0f0")
        v_btn_frame1.pack(pady=2)
        tk.Button(v_btn_frame1, text="Top View", width=12, command=lambda: self.change_view(90, -90)).grid(row=0, column=0, padx=2)
        tk.Button(v_btn_frame1, text="Front View", width=12, command=lambda: self.change_view(0, -90)).grid(row=0, column=1, padx=2)
        tk.Button(v_btn_frame1, text="Side View", width=12, command=lambda: self.change_view(0, 0)).grid(row=0, column=2, padx=2)
        tk.Button(view_frame, text="Default 3D View", width=40, command=lambda: self.change_view(30, -60)).pack(pady=(0, 2))

        ctrl_frame = tk.Frame(view_frame, bg="#f0f0f0")
        ctrl_frame.pack(fill=tk.X, pady=5)
        tk.Button(ctrl_frame, text="▶ Start", font=("Arial", 10, "bold"), bg="#5cb85c", fg="white", command=self.play_anim).pack(side=tk.TOP, fill=tk.X, padx=20, pady=2)
        tk.Button(ctrl_frame, text="⏸ Stop", font=("Arial", 10, "bold"), bg="#f0ad4e", fg="white", command=self.pause_anim).pack(side=tk.TOP, fill=tk.X, padx=20, pady=2)
        tk.Button(ctrl_frame, text="🔄 Reset", font=("Arial", 10, "bold"), bg="#d9534f", fg="white", command=self.reset_anim).pack(side=tk.TOP, fill=tk.X, padx=20, pady=2)

        self.right_frame = tk.Frame(self.root)
        self.right_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True, padx=10, pady=10)
        self.notebook = ttk.Notebook(self.right_frame)
        self.notebook.pack(fill=tk.BOTH, expand=True)
        self.tab_3d = tk.Frame(self.notebook)
        self.tab_sync = tk.Frame(self.notebook)
        self.notebook.add(self.tab_3d, text=" 🚀 3D Trajectory ")
        self.notebook.add(self.tab_sync, text=" 🎯 3D Timing Sync ")

    def populate_tables(self):
        inputs = [
            ("Calculation Mode", self.sd.get('mode', 'N/A').upper()),
            ("Air Pressure", f"{self.sd.get('pressure', 0)} bar"),
            ("Gravity (g)", f"{self.sd.get('g', 9.81):.5f} m/s²"),
            ("Target RPM", f"{self.sd.get('target_rpm', 0)} RPM ({self.sd.get('target_dir', 'CCW')})"),
            ("HW Latency", f"{self.sd.get('hw_latency', 0)} s")
        ]
        calcs = [
            ("Initial Velocity (V0)", f"{self.sd.get('v0', 0):.3f} m/s"),
            ("Time of Flight", f"{self.sd.get('t_flight', 0):.3f} s"),
            ("Firing Delay Needed", f"{self.sd.get('t_delay', 0):.3f} s"),
            ("Calculated Pitch", f"{self.sd.get('pitch', 0):.2f}°"),
            ("Calculated Yaw", f"{self.sd.get('yaw', 0):.2f}°")
        ]
        
        for item in inputs: self.tree_in.insert("", "end", values=item)
        for item in calcs: self.tree_out.insert("", "end", values=item)

        t_flight = self.sd['t_flight']
        t_steps = np.arange(0, t_flight, 0.05)
        if len(t_steps) == 0 or t_steps[-1] < t_flight: t_steps = np.append(t_steps, t_flight)
        
        v0 = self.sd['v0']
        p_rad, y_rad = np.radians(self.sd['pitch']), np.radians(self.sd['yaw'])
        for t in t_steps:
            x = self.sd['p0'][0] + (v0 * np.cos(p_rad) * np.cos(y_rad)) * t
            y = self.sd['p0'][1] + (v0 * np.cos(p_rad) * np.sin(y_rad)) * t
            z = self.sd['p0'][2] + (v0 * np.sin(p_rad)) * t - 0.5 * self.sd['g'] * t**2
            if z < 0 and t == t_flight: z = self.sd['impact_point'][2]
            self.tree_log.insert("", "end", values=(f"{t:.3f}", f"{x:.3f}", f"{y:.3f}", f"{z:.3f}"))

    def style_3d_axes(self, ax):
        ax.zaxis.pane.set_facecolor((0.0, 1.0, 0.0, 0.05)); ax.zaxis.pane.set_edgecolor('green')
        ax.yaxis.pane.set_facecolor((1.0, 0.0, 0.0, 0.05)); ax.yaxis.pane.set_edgecolor('red')
        ax.xaxis.pane.set_facecolor((0.0, 0.0, 1.0, 0.05)); ax.xaxis.pane.set_edgecolor('blue')
        ax.set_xlabel("X (m)", color='blue'); ax.set_ylabel("Y (m)", color='red'); ax.set_zlabel("Z (m)", color='green')
        ax.set_autoscale_on(False) 

    def set_axes_limits(self, ax, all_x, all_y, all_z):
        x_min, x_max = min(all_x), max(all_x)
        y_min, y_max = min(all_y), max(all_y)
        z_min, z_max = min(all_z), max(all_z)
        max_range = max(x_max-x_min, y_max-y_min, z_max-z_min)
        if max_range == 0: max_range = 1.0
        margin = max_range * 0.30
        max_range_with_margin = max_range + (margin*2)
        mid_x, mid_y, mid_z = (x_max + x_min) * 0.5, (y_max + y_min) * 0.5, (z_max + z_min) * 0.5
        half_range = max_range_with_margin * 0.5
        ax.set_xlim(mid_x - half_range, mid_x + half_range)
        ax.set_ylim(mid_y - half_range, mid_y + half_range)
        z_lower = min(0, mid_z - half_range)
        ax.set_zlim(z_lower, z_lower + max_range_with_margin)
        ax.set_box_aspect((1, 1, 1))

    def setup_plots(self):
        t_flight = self.sd['t_flight']
        self.t_arr_3d = np.arange(0, t_flight, 0.05)
        if len(self.t_arr_3d) == 0 or self.t_arr_3d[-1] < t_flight: self.t_arr_3d = np.append(self.t_arr_3d, t_flight)
            
        v0, g = self.sd['v0'], self.sd['g']
        pr, yr = np.radians(self.sd['pitch']), np.radians(self.sd['yaw'])
        p0, impact = self.sd['p0'], self.sd['impact_point']
        
        self.x_3d = p0[0] + (v0 * np.cos(pr) * np.cos(yr)) * self.t_arr_3d
        self.y_3d = p0[1] + (v0 * np.cos(pr) * np.sin(yr)) * self.t_arr_3d
        self.z_3d = p0[2] + (v0 * np.sin(pr)) * self.t_arr_3d - 0.5 * g * self.t_arr_3d**2

        # PRE-CALCULATION FOR SYNC ANIMATION
        self.fps = 20
        self.total_t = self.sd['t_total_sync']
        self.frames_sync = int(np.ceil(self.total_t * self.fps)) + 1
        
        cx, cy, cz = 2.5, 0.0, impact[2] 
        r, math_trigger, math_impact = self.sd['target_radius'], self.sd['math_trigger'], self.sd['math_impact']
        
        self.sync_tx, self.sync_ty, self.sync_tz = np.zeros(self.frames_sync), np.zeros(self.frames_sync), np.full(self.frames_sync, cz)
        self.sync_bx, self.sync_by, self.sync_bz = np.zeros(self.frames_sync), np.zeros(self.frames_sync), np.zeros(self.frames_sync)
        self.sync_b_vis, self.sync_t_blink = np.zeros(self.frames_sync, dtype=bool), np.zeros(self.frames_sync, dtype=bool)
        self.sync_titles, self.sync_colors = [], []
        
        t_delay, t_hw = self.sd['t_delay'], self.sd['hw_latency']
        t_launch = t_delay + t_hw
        
        time_per_rev = 360.0 / abs(self.sd['omega'])

        for frame in range(self.frames_sync):
            raw_t = frame / self.fps
            
            # --- FIX ERROR 7.5 DEGREE --- 
            # ล็อคเวลาเฟรมสุดท้ายให้เท่ากับ total_t พอดีเป๊ะ ไม่ให้เกินเด็ดขาด
            if raw_t >= self.total_t:
                t = self.total_t
                is_impact_frame = True
            else:
                t = raw_t
                is_impact_frame = False
            
            curr_ang = math_trigger + (self.sd['omega'] * t)
            self.sync_tx[frame] = cx + r * math.cos(math.radians(curr_ang))
            self.sync_ty[frame] = cy + r * math.sin(math.radians(curr_ang))
            
            time_since_trigger = t % time_per_rev
            self.sync_t_blink[frame] = True if time_since_trigger < 0.15 else False
            
            if is_impact_frame:
                self.sync_b_vis[frame] = True
                # ล็อกตำแหน่งบอลให้เป๊ะที่จุดตก
                self.sync_bx[frame] = p0[0] + (v0 * math.cos(pr) * math.cos(yr)) * t_flight
                self.sync_by[frame] = p0[1] + (v0 * math.cos(pr) * math.sin(yr)) * t_flight
                self.sync_bz[frame] = cz
                
                # Validation จะใช้มุมที่เวลาชนแบบ Exact Math! ไม่มีเศษ Error แน่นอน
                msg, color, err = TimingLogic.validate_impact(curr_ang, math_impact, self.sd['tolerance'])
                self.sync_titles.append(f"Status: IMPACT! [{msg}]\nError: {err:.1f} deg")
                self.sync_colors.append(color)
            elif t < t_delay:
                self.sync_b_vis[frame] = False
                self.sync_titles.append(f"Status: DELAYING (Waiting to fire)\nTime: {t:.2f} s")
                self.sync_colors.append("orange")
            elif t < t_launch:
                self.sync_b_vis[frame] = True
                self.sync_bx[frame], self.sync_by[frame], self.sync_bz[frame] = p0[0], p0[1], p0[2]
                self.sync_titles.append(f"Status: HW FIRING (Valve Opening)\nTime: {t:.2f} s")
                self.sync_colors.append("red")
            else:
                self.sync_b_vis[frame] = True
                t_proj = min(t - t_launch, t_flight)
                self.sync_bx[frame] = p0[0] + (v0 * math.cos(pr) * math.cos(yr)) * t_proj
                self.sync_by[frame] = p0[1] + (v0 * math.cos(pr) * math.sin(yr)) * t_proj
                self.sync_bz[frame] = p0[2] + (v0 * math.sin(pr)) * t_proj - 0.5 * g * t_proj**2
                if self.sync_bz[frame] < 0 and t_proj >= t_flight: self.sync_bz[frame] = cz
                self.sync_titles.append(f"Status: PROJECTILE IN FLIGHT!\nTime: {t:.2f} s")
                self.sync_colors.append("blue")

        # TAB 1: 3D
        self.fig_3d = plt.figure(figsize=(6, 5))
        self.ax_3d = self.fig_3d.add_subplot(111, projection='3d')
        self.canvas_3d = FigureCanvasTkAgg(self.fig_3d, master=self.tab_3d)
        self.canvas_3d.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        self.style_3d_axes(self.ax_3d)
        
        self.line_3d, = self.ax_3d.plot([], [], [], lw=2.5, color='blue', label='Trajectory')
        self.point_3d, = self.ax_3d.plot([], [], [], 'ro', markersize=8)
        self.ax_3d.scatter(*p0, color='green', s=50, label='Start Point')
        self.ax_3d.scatter(*impact, color='red', marker='x', s=100, label='Impact Point')
        self.ax_3d.legend()
        
        # TAB 2: Sync
        self.fig_sync = plt.figure(figsize=(6, 5))
        self.ax_sync = self.fig_sync.add_subplot(111, projection='3d')
        self.canvas_sync = FigureCanvasTkAgg(self.fig_sync, master=self.tab_sync)
        self.canvas_sync.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        self.style_3d_axes(self.ax_sync)
        
        self.status_text = self.ax_sync.text2D(0.02, 0.95, "", transform=self.ax_sync.transAxes, fontsize=11, weight='bold', bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))
        
        theta = np.linspace(0, 2*np.pi, 100)
        circle_x, circle_y, circle_z = cx + r * np.cos(theta), cy + r * np.sin(theta), np.full_like(theta, cz)
        self.ax_sync.plot(circle_x, circle_y, circle_z, color='gray', linestyle='--', label='Target Path')
        self.ax_sync.plot([cx], [cy], [cz], 'k^', markersize=10, label='Motor Center (2.5, 0)')
        
        trig_rad_plot = math.radians(math_trigger)
        self.ax_sync.plot([cx, cx + r * math.cos(trig_rad_plot)], [cy, cy + r * math.sin(trig_rad_plot)], [cz, cz], color='orange', linestyle=':', label='Trigger Line')
        
        self.ax_sync.scatter(*p0, color='green', s=50, label='Muzzle')
        self.ax_sync.scatter(*impact, color='red', marker='x', s=100, label='Impact Point')
        
        self.target_sync, = self.ax_sync.plot([], [], [], 'bo', markersize=10, label='Rotating Target')
        self.ball_sync, = self.ax_sync.plot([], [], [], 'go', markersize=8, label='Projectile')
        self.ax_sync.legend()

        all_x_lim = np.concatenate([self.x_3d, circle_x, [cx]])
        all_y_lim = np.concatenate([self.y_3d, circle_y, [cy]])
        all_z_lim = np.concatenate([self.z_3d, circle_z, [cz]])
        self.set_axes_limits(self.ax_3d, all_x_lim, all_y_lim, all_z_lim)
        self.set_axes_limits(self.ax_sync, all_x_lim, all_y_lim, all_z_lim)
        self.change_view(30, -60) 
        
        self.anim_3d = FuncAnimation(self.fig_3d, self.update_3d, frames=len(self.t_arr_3d), interval=50, blit=False, repeat=True)
        self.anim_sync = FuncAnimation(self.fig_sync, self.update_sync, frames=self.frames_sync, interval=50, blit=False, repeat=True)
        self.pause_anim()

    def update_monitor(self, t_global, t_proj, state):
        self.lbl_mon_time.config(text=f"Time: {t_global:.3f} s | State: {state}")
        
        v0, g = self.sd['v0'], self.sd['g']
        pr, yr = np.radians(self.sd['pitch']), np.radians(self.sd['yaw'])
        p0 = self.sd['p0']
        
        vx_c = v0 * math.cos(pr) * math.cos(yr)
        vy_c = v0 * math.cos(pr) * math.sin(yr)
        vz0 = v0 * math.sin(pr)
        
        if state in ["DELAY", "HW LATENCY"]:
            x, y, z = p0
            vx, vy, vz, ax, ay, az = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
        elif state == "FLYING":
            x = p0[0] + vx_c * t_proj
            y = p0[1] + vy_c * t_proj
            z = p0[2] + vz0 * t_proj - 0.5 * g * t_proj**2
            vx, vy, vz = vx_c, vy_c, vz0 - g * t_proj
            ax, ay, az = 0.0, 0.0, -g
        else:
            x = p0[0] + vx_c * t_proj
            y = p0[1] + vy_c * t_proj
            z = p0[2] + vz0 * t_proj - 0.5 * g * t_proj**2
            if z < self.sd['impact_point'][2]: z = self.sd['impact_point'][2]
            vx, vy, vz, ax, ay, az = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
            
        self.tree_mon.item("X", values=("X", f"{x:.3f}", f"{vx:.3f}", f"{ax:.3f}"))
        self.tree_mon.item("Y", values=("Y", f"{y:.3f}", f"{vy:.3f}", f"{ay:.3f}"))
        self.tree_mon.item("Z", values=("Z", f"{z:.3f}", f"{vz:.3f}", f"{az:.3f}"))

    def update_3d(self, i):
        t = self.t_arr_3d[i]
        self.line_3d.set_data(self.x_3d[:i+1], self.y_3d[:i+1])
        self.line_3d.set_3d_properties(self.z_3d[:i+1])
        self.point_3d.set_data([self.x_3d[i]], [self.y_3d[i]])
        self.point_3d.set_3d_properties([self.z_3d[i]])
        
        if self.notebook.index(self.notebook.select()) == 0:
            state = "FLYING" if t < self.sd['t_flight'] else "IMPACT"
            self.update_monitor(t, t, state)
            
        return self.line_3d, self.point_3d

    def update_sync(self, frame):
        raw_t = frame / self.fps
        t = self.total_t if raw_t >= self.total_t else raw_t
        
        self.target_sync.set_data([self.sync_tx[frame]], [self.sync_ty[frame]])
        self.target_sync.set_3d_properties([self.sync_tz[frame]])
        
        if self.sync_t_blink[frame]:
            self.target_sync.set_color('#FFD700')
            self.target_sync.set_markersize(15)
        else:
            self.target_sync.set_color('blue')
            self.target_sync.set_markersize(10)
        
        if self.sync_b_vis[frame]:
            self.ball_sync.set_data([self.sync_bx[frame]], [self.sync_by[frame]])
            self.ball_sync.set_3d_properties([self.sync_bz[frame]])
        else:
            self.ball_sync.set_data([], [])
            self.ball_sync.set_3d_properties([])
            
        title = self.sync_titles[frame]
        color = self.sync_colors[frame]
        self.status_text.set_text(title)
        self.status_text.set_color(color)
        
        if self.notebook.index(self.notebook.select()) == 1:
            t_launch = self.sd['t_delay'] + self.sd['hw_latency']
            t_flight = self.sd['t_flight']
            
            if t < self.sd['t_delay']: state = "DELAY"
            elif t < t_launch: state = "HW LATENCY"
            elif t < self.total_t: state = "FLYING" 
            else: state = "IMPACT"
            
            t_proj = max(0, min(t - t_launch, t_flight))
            if state in ["DELAY", "HW LATENCY"]: t_proj = 0
            
            self.update_monitor(t, t_proj, state)
        
        if 'IMPACT' in title:
            self.anim_sync.event_source.stop()
            
        return self.target_sync, self.ball_sync

    def change_view(self, elev, azim):
        self.ax_3d.view_init(elev=elev, azim=azim)
        self.ax_sync.view_init(elev=elev, azim=azim)
        self.canvas_3d.draw()
        self.canvas_sync.draw()

    def play_anim(self):
        if not self.is_playing:
            if self.anim_3d and self.anim_3d.event_source: self.anim_3d.resume()
            if self.anim_sync and self.anim_sync.event_source: self.anim_sync.resume()
            self.is_playing = True

    def pause_anim(self):
        if self.is_playing:
            if self.anim_3d and self.anim_3d.event_source: self.anim_3d.pause()
            if self.anim_sync and self.anim_sync.event_source: self.anim_sync.pause()
            self.is_playing = False

    def reset_anim(self):
        self.pause_anim()
        if self.anim_3d and self.anim_3d.event_source: self.anim_3d.event_source.stop()
        if self.anim_sync and self.anim_sync.event_source: self.anim_sync.event_source.stop()
        
        self.anim_3d = FuncAnimation(self.fig_3d, self.update_3d, frames=len(self.t_arr_3d), interval=50, blit=False, repeat=True)
        self.anim_sync = FuncAnimation(self.fig_sync, self.update_sync, frames=self.frames_sync, interval=50, blit=False, repeat=True)
        
        self.update_3d(0); self.canvas_3d.draw()
        self.update_sync(0); self.canvas_sync.draw()
        self.pause_anim()

    def run(self): self.root.mainloop()

if __name__ == "__main__":
    app = SimulationWindow()
    app.run()