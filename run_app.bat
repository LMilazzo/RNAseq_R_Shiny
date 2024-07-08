@echo off

REM pull latest changes
git pull origin main

REM run shiny
Rscript run_app.R

pause