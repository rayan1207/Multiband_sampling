import random
import matplotlib.pyplot as plt
import numpy as np




# def fun(x,y):
#     return np.cos(x)-np.cos(y)


# print(fun(-2.42807 ,-0.390257))
p = np.pi
vertices = [(-2*p/3, 2*p/np.sqrt(3)), (2*p/3, 2*p/np.sqrt(3)), (4*p/3, 0), (2*p/3, -2*p/np.sqrt(3)), (-2*p/3, -2*p/np.sqrt(3)), (-4*p/3, 0)]

# Function to check if a point is inside the hexagon using the Ray Casting method
def is_inside(point):
    x, y = point
    odd_nodes = False
    j = len(vertices) - 1

    for i in range(len(vertices)):
        xi, yi = vertices[i]
        xj, yj = vertices[j]
        if yi < y and yj >= y or yj < y and yi >= y:
            if xi + (y - yi) / (yj - yi) * (xj - xi) < x:
                odd_nodes = not odd_nodes
        j = i

    return odd_nodes

# Find the bounding box of the hexagon
min_x = min(x for x, y in vertices)
max_x = max(x for x, y in vertices)
min_y = min(y for x, y in vertices)
max_y = max(y for x, y in vertices)

# Generate random points within the bounding box
total_points = 100000# You can increase this number for better accuracy
points_inside_hexagon = 0

# Lists to store the x and y coordinates of the points inside and outside the hexagon
x_inside, y_inside = [], []
x_outside, y_outside = [], []

for _ in range(total_points):
    x_rand = random.uniform(min_x, max_x)
    y_rand = random.uniform(min_y, max_y)
    if is_inside((x_rand, y_rand)):
        points_inside_hexagon += 1
        x_inside.append(x_rand)
        y_inside.append(y_rand)
    else:
        x_outside.append(x_rand)
        y_outside.append(y_rand)

# Calculate the area using the Monte Carlo method
area_bounding_box = (max_x - min_x) * (max_y - min_y)
area_hexagon = (points_inside_hexagon / total_points) * area_bounding_box

print("Estimated area of the hexagon:", area_hexagon)

# Plot the vertices and random points
plt.figure(figsize=(8, 6))
plt.scatter(*zip(*vertices), color='red', label='Vertices')
plt.scatter(x_inside, y_inside, color='green', alpha=0.5, label='Points Inside Hexagon')
plt.scatter(x_outside, y_outside, color='blue', alpha=0.5, label='Points Outside Hexagon')
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Scatter Plot of Hexagon')
plt.legend()
plt.grid(True)
plt.show()
