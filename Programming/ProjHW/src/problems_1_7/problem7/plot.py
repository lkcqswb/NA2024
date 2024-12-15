import sys
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import os
import json

def read_config(file_name):
    try:
        with open(file_name, 'r') as file:
            config = json.load(file)
            return config
    except Exception as e:
        print(f"Error reading file: {e}")
        return None

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


def plot_sphere(ax, center, radius):
    # 生成球面坐标
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)
    x = center[0] + radius * np.outer(np.cos(u), np.sin(v))
    y = center[1] + radius * np.outer(np.sin(u), np.sin(v))
    z = center[2] + radius * np.outer(np.ones(np.size(u)), np.cos(v))
    
    # 绘制球面
    ax.plot_surface(x, y, z, color='b', alpha=0.1)  # alpha 控制透明度

def plot_points(points, title, center, radius):
    if not points:
        print("No points to plot.")
        return

    # 将点分为 x, y 和 z 坐标
    x, y, z = zip(*points)

    # 创建三维坐标轴
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # 绘制在球面上的细线
    ax.plot(x, y, z, color='red', linewidth=1)  # 设置 linewidth 为 1

    ax.set_title(title)
    ax.set_xlabel('X-axis')
    ax.set_ylabel('Y-axis')
    ax.set_zlabel('Z-axis')

    # 绘制球面
    plot_sphere(ax, center, radius)

    # 显示图形
    plt.savefig(title+".png")

if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Usage: python plot.py <filename> <radius> <center_x> <center_y> <center_z>")
        sys.exit(1)

    file_name = sys.argv[1]
    radius = float(sys.argv[2])
    center = [float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5])]
    
    
    config = read_config("config.json")
    if config is None:
        sys.exit(1)
    center = config["centre"]
    radius = config["radius"]
    points = config["points"]
    title = "3D Sphere with Curve"
    plot_points(points, title, center, radius)

    points = read_points(file_name)
    title = os.path.splitext(file_name)[0]
    plot_points(points, title, center, radius)