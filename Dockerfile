# Use an official Python runtime as a parent image
FROM python:3.14.0a3-slim

# Set the working directory in the container
WORKDIR /app

# Copy the current directory contents into the container at /app

COPY . /app

# Install any necessary system dependencies (optional, adjust as needed)
# RUN apt-get update && apt-get install -y <your-system-dependencies>

# Install the Python package using setup.py
RUN pip install .

# Expose any necessary ports (if your application needs it, e.g., for Dash)
EXPOSE 8050

USER 1001

# Run the application (replace this with your application command)
ENV HOST=0.0.0.0
CMD ["gunicorn", "arcs.dash_app:server", "-b", ":8050"]
