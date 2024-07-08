#!/bin/bash

# Pull the latest changes from the GitHub repository
git pull origin main

# Run the Shiny app
Rscript run_app.R