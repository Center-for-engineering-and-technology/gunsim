import pyautogui
import time

print("Координаты мыши (нажмите Ctrl+C для остановки):")
try:
    while True:
        x, y = pyautogui.position()
        print(f"X: {x:4} | Y: {y:4}", end='\r')
        time.sleep(0.1)
except KeyboardInterrupt:
    print("\nГотово!")