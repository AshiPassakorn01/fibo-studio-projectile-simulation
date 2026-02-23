import numpy as np
import math

class PhysicsEngine:
    @staticmethod
    def calculate_v0(p_bar, bore_mm, stroke_mm, mass_kg):
        p_pa = p_bar * 100000
        r_m = (bore_mm / 2) / 1000
        area = np.pi * (r_m ** 2)
        stroke_m = stroke_mm / 1000
        force = p_pa * area
        work = force * stroke_m
        v0 = np.sqrt((2 * work) / mass_kg)
        return v0, work, force

    @staticmethod
    def forward_kinematics(v0, pitch_deg, yaw_deg, p0, target_z, g):
        pitch, yaw = np.radians(pitch_deg), np.radians(yaw_deg)
        vx, vy, vz = v0 * np.cos(pitch) * np.cos(yaw), v0 * np.cos(pitch) * np.sin(yaw), v0 * np.sin(pitch)
        a, b, c = 0.5 * g, -vz, target_z - p0[2]
        
        discriminant = b**2 - 4*a*c
        if discriminant < 0: raise ValueError("Projectile cannot reach the specified Z level.")
            
        t_flight = max((-b + np.sqrt(discriminant)) / (2*a), (-b - np.sqrt(discriminant)) / (2*a))
        if t_flight <= 0: raise ValueError("Invalid target Z calculation.")
            
        impact_x, impact_y = p0[0] + vx * t_flight, p0[1] + vy * t_flight
        distance = np.sqrt((impact_x - p0[0])**2 + (impact_y - p0[1])**2)
        
        return {'impact_point': (impact_x, impact_y, target_z), 't_flight': t_flight, 'distance': distance}

    @staticmethod
    def inverse_kinematics(v0, p0, p_target, g):
        dx, dy, dz = p_target[0] - p0[0], p_target[1] - p0[1], p_target[2] - p0[2]
        yaw = np.degrees(np.arctan2(dy, dx))
        dist_h = np.sqrt(dx**2 + dy**2)
        
        A = (g * dist_h**2) / (2 * v0**2)
        discriminant = dist_h**2 - 4 * A * (A + dz)
        
        if discriminant < 0: raise ValueError("Target is too far for the given air pressure.")
            
        u1, u2 = (dist_h + np.sqrt(discriminant)) / (2 * A), (dist_h - np.sqrt(discriminant)) / (2 * A)
        valid_pitches = [p for p in (np.degrees(np.arctan(u1)), np.degrees(np.arctan(u2))) if 40 < p < 90]
        if not valid_pitches: raise ValueError("No valid pitch angle found (must be between 40-90 degrees).")
            
        best_pitch = max(valid_pitches)
        t_flight = dist_h / (v0 * np.cos(np.radians(best_pitch)))
        
        return {'pitch': best_pitch, 'yaw': yaw, 't_flight': t_flight, 'distance': dist_h}

class TimingLogic:
    @staticmethod
    def calculate_delay(rpm, trigger_angle, impact_angle, t_hw, t_flight, direction="CCW"):
        if rpm <= 0: raise ValueError("RPM must be greater than 0")

        omega_deg_per_sec = rpm * 6.0 
        time_per_rev = 360.0 / omega_deg_per_sec
        
        # ค้นหาระยะทางเชิงมุมตามทิศทางหมุน
        if direction.upper() == "CW":
            delta_theta = (trigger_angle - impact_angle) % 360
            omega_anim = -omega_deg_per_sec
        else: # CCW
            delta_theta = (impact_angle - trigger_angle) % 360
            omega_anim = omega_deg_per_sec
            
        t_target_base = delta_theta / omega_deg_per_sec
        
        # เงื่อนไข: หากจุด Trigger อยู่ใกล้เป้าเกินไปจนใช้เวลาน้อยกว่า 1 วินาที ให้รอรอบหน้าเลย
        if t_target_base < 1.0:
            t_target_base += time_per_rev
            
        total_projectile_time = t_hw + t_flight
        
        if total_projectile_time <= t_target_base:
            n = 0
        else:
            n = math.ceil((total_projectile_time - t_target_base) / time_per_rev)
            
        t_target_actual = t_target_base + (n * time_per_rev)
        t_delay = t_target_actual - total_projectile_time
        
        return t_delay, t_target_actual, omega_anim

    @staticmethod
    def validate_impact(target_angle_at_impact, impact_angle, tolerance=5.0):
        error = (target_angle_at_impact - impact_angle) % 360
        if error > 180: error -= 360
            
        if abs(error) <= tolerance: return "HIT", "green", error
        elif error > 0: return "LATE", "red", error
        else: return "EARLY", "orange", error