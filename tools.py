import json
import os

DATA_FILE = "data.json"

def load_data():
    """โหลดข้อมูลจากไฟล์ JSON"""
    if not os.path.exists(DATA_FILE):
        return {}
    with open(DATA_FILE, 'r', encoding='utf-8') as f:
        try:
            return json.load(f)
        except json.JSONDecodeError:
            return {}

def update_data(new_data):
    """อัปเดตข้อมูลใหม่ลงในไฟล์ JSON"""
    data = load_data()
    data.update(new_data)
    with open(DATA_FILE, 'w', encoding='utf-8') as f:
        json.dump(data, f, indent=4)