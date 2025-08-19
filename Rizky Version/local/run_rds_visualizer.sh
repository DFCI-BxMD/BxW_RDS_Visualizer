#!/bin/sh

echo "Starting Single-Cell Analyzer..."

docker load -i rds_vis_local.tar

# Run the Docker container
docker run --rm -d \
  -p 3838:3838 \
  -v $HOME:/srv/shiny-server/data \
  --name single_cell_analyzer \
  rds_vis_local:latest

sleep 15

# Detect the operating system and use the correct command
echo "Opening browser to http://localhost:3838..."

if [[ "$OSTYPE" == "linux-gnu"* ]]; then
    xdg-open http://localhost:3838
elif [[ "$OSTYPE" == "darwin"* ]]; then
    open http://localhost:3838
elif [[ "$OSTYPE" == "cygwin" || "$OSTYPE" == "msys" ]]; then
    start http://localhost:3838
else
    echo "Could not detect OS. Please open http://localhost:3838 manually."
fi

echo "App is running. To stop it, run: docker stop single_cell_analyzer"
echo "The container will automatically stop after 3 hours."

(
sleep 10800
docker stop rds_visualizer
) &

echo "Container has been stopped. Please run the script again to start a new session"