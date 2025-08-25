import numpy as np
import math

# -------------------------------
# Tool Funcitons
# -------------------------------
def wrap_pi(a):
    """把角度限制在 (-pi, pi]"""
    return (a + np.pi) % (2*np.pi) - np.pi

def azimuth(x1, y1, x2, y2):
    """計算從 (x1,y1) 指向 (x2,y2) 的方位角 (弧度)"""
    return math.atan2(y2 - y1, x2 - x1)

# -------------------------------
# 四個已知點
# -------------------------------
known = {
    "A": (0.0, 0.0),
    "B": (8.0, 0.6),
    "D": (7.2, 6.2),
    "E": (1.2, 7.4),
}

# -------------------------------
# 假設一個真實未知點
# -------------------------------
X_true, Y_true = 4.3, 3.9

# 定義每個站的參考線 (自己決定)
Ref = {"A": "B", "B": "A", "D": "B", "E": "A"}

# -------------------------------
# 生成模擬觀測
# -------------------------------
obs = []
for st, (Xi, Yi) in known.items():
    # 距離
    dist = math.hypot(X_true - Xi, Y_true - Yi)
    # 方位角差 (Az(未知) - Az(參考線))
    refX, refY = known[Ref[st]]
    az_ref = azimuth(Xi, Yi, refX, refY)
    az_tar = azimuth(Xi, Yi, X_true, Y_true)
    ang = wrap_pi(az_tar - az_ref)
    obs.append({"station": st, "dist": dist, "angle": ang})

# -------------------------------
# 最小二乘 假設初始值 4.1,3.4
# -------------------------------
def solve_unknown(known, obs, Ref, X0=4.1, Y0=3.4, max_iter=20, tol=1e-10):
    X, Y = X0, Y0
    for it in range(max_iter):
        w_list = []
        A_list = []
        for o in obs:
            st = o["station"]
            Xi, Yi = known[st]
            dX, dY = X - Xi, Y - Yi
            r = math.hypot(dX, dY)

            # 距離方程
            # w = (obs - calc)
            w_d = o["dist"] - r
            # 偏導
            a_dX = dX / r
            a_dY = dY / r
            A_list.append([a_dX, a_dY])
            w_list.append(w_d)

            # 角度方程
            phi = math.atan2(dY, dX)
            refX, refY = known[Ref[st]]
            Theta = azimuth(Xi, Yi, refX, refY)
            w_th = wrap_pi(o["angle"] - (phi - Theta))
            # 偏導
            a_thX = dY / (r**2)
            a_thY = -dX / (r**2)
            A_list.append([a_thX, a_thY])
            w_list.append(w_th)

        A = np.array(A_list)
        w = np.array(w_list).reshape(-1,1)

        # 正規方程
        N = A.T @ A
        u = A.T @ w
        dx = np.linalg.solve(N, u)

        X += dx[0,0]
        Y += dx[1,0]

        if np.linalg.norm(dx) < tol:
            break

    v = A @ dx - w   # 殘差
    dof = len(w) - 2
    sigma0 = math.sqrt(float(v.T @ v) / dof)
    Cov = sigma0**2 * np.linalg.inv(N)
    sigX = math.sqrt(Cov[0,0])
    sigY = math.sqrt(Cov[1,1])
    rho = Cov[0,1]/(sigX*sigY)

    return {
        "X": X, "Y": Y,
        "sigma0": sigma0,
        "sigma_X": sigX, "sigma_Y": sigY, "rho": rho,
        "iterations": it+1
    }

# -------------------------------
# 執行
# -------------------------------
ans = solve_unknown(known, obs, Ref, X0=2, Y0=2)
print("解算結果：", ans)
print("真實值： X=%.3f, Y=%.3f" % (X_true, Y_true))
