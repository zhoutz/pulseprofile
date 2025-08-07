import numpy as np
import matplotlib.pyplot as plt

# 创建3D图形
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection="3d")

# 生成球面网格
u = np.linspace(0, 2 * np.pi, 100)  # 经度角
v = np.linspace(0, np.pi, 50)  # 纬度角
x = np.outer(np.cos(u), np.sin(v))  # 球面x坐标
y = np.outer(np.sin(u), np.sin(v))  # 球面y坐标
z = np.outer(np.ones(np.size(u)), np.cos(v))  # 球面z坐标

# 绘制半透明球面
ax.plot_surface(x, y, z, color="b", alpha=0.1)

# 绘制经度线（固定纬度，改变经度）
for i in range(0, 100, 10):
    ax.plot(x[i, :], y[i, :], z[i, :], color="black", linewidth=0.5)

# 绘制纬度线（固定经度，改变纬度）
for j in range(0, 50, 5):
    ax.plot(x[:, j], y[:, j], z[:, j], color="black", linewidth=0.5)

# 示例点集（球坐标：(theta, phi)）
# points = [
#     (np.pi/4, np.pi/3),     # 点1
#     (np.pi/3, np.pi/2),     # 点2
#     (np.pi/2, 2*np.pi/3),   # 点3
#     (2*np.pi/3, 5*np.pi/4), # 点4
#     (3*np.pi/4, 3*np.pi/2)  # 点5
# ]
points = np.loadtxt("grid.txt")
print(points.shape)
points[:, 0] = np.arccos(points[:, 0])


# 将球坐标转换为笛卡尔坐标
def spherical_to_cartesian(theta, phi):
    x = np.sin(theta) * np.cos(phi)
    y = np.sin(theta) * np.sin(phi)
    z = np.cos(theta)
    return x, y, z


# 绘制点
xs, ys, zs = [], [], []
for theta, phi in points:
    x, y, z = spherical_to_cartesian(theta, phi)
    xs.append(x)
    ys.append(y)
    zs.append(z)

ax.scatter(xs, ys, zs, color="red", s=1, depthshade=True)

# 添加标签和标题
ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")
ax.set_title("3D Unit Sphere with Points")
ax.set_box_aspect([1, 1, 1])  # 等比例缩放坐标轴

plt.show()
