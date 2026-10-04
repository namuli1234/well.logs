@echo off
REM Builds WellLogPro.exe (requires Python 3.9+ from python.org, "Add to PATH" ticked)
python -m pip install --upgrade pip
python -m pip install numpy matplotlib pyinstaller
python -m PyInstaller --noconfirm --onefile --windowed --name WellLogPro app.py
echo.
echo Done. Your app is: dist\WellLogPro.exe
pause
