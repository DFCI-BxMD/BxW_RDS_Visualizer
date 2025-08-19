#!/bin/sh

Write-Host "Starting Single-Cell Analyzer..."

# Run the Docker container
docker run --rm -d `
  -p 3838:3838 `
  -v $HOME:/srv/shiny-server/data `
  --name single_cell_analyzer `
  rds_vis_local:latest

Start-Sleep -Seconds 15

Write-Host "Opening browser to http://localhost:3838..."

Start-Process http://localhost:3838

Write-Host "App is running. To stop it, run: docker stop single_cell_analyzer"
Write-Host "The container will automatically stop after 3 hours."

Start-Sleep -Seconds 10800
docker stop rds_visualizer

Write-Host "Container has been stopped. Please run the script again to start a new session"