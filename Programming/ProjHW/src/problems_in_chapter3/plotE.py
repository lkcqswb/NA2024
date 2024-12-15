import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def parse_data(filename):
    with open(filename, 'r', encoding='utf-8') as f:
        lines = f.readlines()
    
    figures = {'fig1': {}, 'fig2': {}, 'fig3': {}}
    current_figure = None
    

    for line in lines:
        line = line.strip()
        if not line:
            continue
        
        # 处理图形标签
        if line.endswith('fig1'):
            current_figure = 'fig1'

            continue
        elif line.endswith('fig2'):
            current_figure = 'fig2'

            continue
        elif line.endswith('fig3'):
            current_figure = 'fig3'
            continue
        
        # 如果还没有确定图形,跳过该行
        if current_figure is None:
            continue
        
        # Parse curve data
        if line.startswith('equal length ppForm'):
            curve_name = 'equal length ppForm'
            parts = line.split(',')
            coords = list(map(float, parts[1:]))
            figures[current_figure][curve_name] = figures[current_figure].get(curve_name, []) + [coords]
            
        elif line.startswith('equal length BSpline'):
            curve_name = 'equal length BSpline'
            parts = line.split(',')
            coords = list(map(float, parts[1:]))
            figures[current_figure][curve_name] = figures[current_figure].get(curve_name, []) + [coords]
            
        elif line.startswith('cumulative chordal length ppform'):
            curve_name = 'cumulative chordal length ppform'
            parts = line.split(',')
            coords = list(map(float, parts[1:]))
            figures[current_figure][curve_name] = figures[current_figure].get(curve_name, []) + [coords]
            
        elif line.startswith('cumulative chordal length BSpline'):
            curve_name = 'cumulative chordal length BSpline'
            parts = line.split(',')
            coords = list(map(float, parts[1:]))
            figures[current_figure][curve_name] = figures[current_figure].get(curve_name, []) + [coords]
            
        else:
            print(f"忽略行: {line}")
    

    
    return figures

def plot_figures(data,file_name):
    plt.figure(figsize=(16, 5))
    
    # Figure 1: 2D plot
    plt.subplot(131)
    for curve_name, points in data['fig1'].items():
        x = [p[0] for p in points]
        y = [p[1] for p in points]
        plt.plot(x, y, label=curve_name, marker='o')
    plt.title('Figure 1 (2D)')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.legend()
    plt.grid(True)

    # Figure 2: 2D plot
    plt.subplot(132)
    for curve_name, points in data['fig2'].items():
        x = [p[0] for p in points]
        y = [p[1] for p in points]
        plt.plot(x, y, label=curve_name, marker='o')
    plt.title('Figure 2 (2D)')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.legend()
    plt.grid(True)

    # Figure 3: 3D plot
    ax = plt.subplot(133, projection='3d')
    for curve_name, points in data['fig3'].items():
        x = [p[0] for p in points]
        y = [p[1] for p in points]
        z = [p[2] for p in points]
        ax.plot(x, y, z, label=curve_name, marker='o')
    ax.set_title('Figure 3 (3D)')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.legend()

    plt.tight_layout()
    plt.savefig(file_name+".png")

import sys
if len(sys.argv) != 2:
        print("Usage: python plot_points.py <filename>")
        sys.exit(1)

file_name = sys.argv[1]
data = parse_data('Pe.txt')
plot_figures(data,file_name)