import pandas as pd
import matplotlib.pyplot as plt
import math

def calculate_distance(point1, point2):
    return math.sqrt((point1[0] - point2[0])**2 + (point1[1] - point2[1])**2)

def calculate_connections(points, threshold_radius):
    connections = []
    for point in points:
        count = sum(calculate_distance(point, other_point) < threshold_radius for other_point in points if point != other_point)
        connections.append(count)
    return connections

excel_file_path = 'Champ de pissenlits et de sauge des pres.xlsx'
df = pd.read_excel(excel_file_path)
flowers = list(zip(df['x'], df['y']))
hive_location = (500, 500)
threshold_radius = 100  # Threshold for circular area

# Calculate the number of connections for each flower
connections = calculate_connections(flowers, threshold_radius)

# Normalize connections for color mapping
max_connections = max(connections)
normalized_connections = [c/max_connections for c in connections]

# Plotting
plt.figure(figsize=(12, 10))
x_coords, y_coords = zip(*flowers)  # Extract x and y coordinates

# Draw flowers with color based on number of connections
for i, point in enumerate(flowers):
    plt.scatter(*point, color=plt.cm.Reds(normalized_connections[i]))

plt.scatter(*hive_location, color='red', label='Hive', marker='x')

# Drawing circles around each point based on the threshold
for point in flowers:
    circle = plt.Circle(point, threshold_radius, color='gray', fill=False, linestyle='--', linewidth=1)
    plt.gca().add_patch(circle)

plt.title('Flowers with Connectivity Visualization')
plt.xlabel('X Coordinate')
plt.ylabel('Y Coordinate')
plt.legend()
plt.grid(True)
plt.show()
