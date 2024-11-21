import sys
import matplotlib.pyplot as plt
import os
from mpl_toolkits.mplot3d import Axes3D 

def read_points(file_name):
    points = []
    try:
        with open(file_name, 'r') as file:
            for line in file:
                point = list(map(float, line.strip().split()))
                points.append(point)
    except Exception as e:
        print(f"Error reading file: {e}")
    return points

def plot_lines(points, title, dimensions):
    if not points:
        print("No points to plot.")
        return

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d') 

    x = [point[0] for point in points]
    
    for i in range(1, dimensions + 1):
        y = [point[i] for point in points]
        ax.plot(x, y, zs=i, zdir='y', label=f'Dimension {i}')

    ax.set_title(title)
    ax.set_xlabel('t-axis')
    ax.set_ylabel('Dimensions')
    ax.set_zlabel('Y-axis')
    ax.legend() 
    plt.show()

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python plot_points.py <filename> <dimensions>")
        sys.exit(1)

    file_name = sys.argv[1]
    dimensions = int(sys.argv[2])  
    points = read_points(file_name)
    title = os.path.splitext(file_name)[0]
    plot_lines(points, title, dimensions)