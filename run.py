import numpy as np
import math

# -------------------------------
# 工具函式
# -------------------------------
def wrap_pi(a):
    """把角度限制在 (-pi, pi]"""
    return (a + np.pi) % (2*np.pi) - np.pi

def azimuth(x1, y1, x2, y2):
    """方位角 (弧度)"""
    return math.atan2(y2 - y1, x2 - x1)

def deg2rad(x_deg):
    return x_deg * math.pi / 180.0

# -------------------------------
# 四個已知點
# -------------------------------
known = {
    "A": (0.0, 0.0),
    "B": (8.0, 0.6),
    "D": (7.2, 6.2),
    "E": (1.2, 7.4),
}

# 每站的參考線終點（量角的起始方向 i->Ref[i]）
# 例如在 A 站用 AB 為參考線 → Ref["A"]="B"
Ref = {"A": "B", "B": "A", "D": "B", "E": "A"}

# -------------------------------
# 將「使用者輸入的距離與角度」轉成 obs
# distances: dict，如 {"A": sA, "B": sB, "D": sD, "E": sE} （單位：公尺）
# angles:    dict，如 {"A": aA, "B": aB, "D": aD, "E": aE} （單位：角度或弧度；由 angle_unit 指定）
# angle_unit: "deg" 或 "rad"
# -------------------------------
def build_obs_from_user(known, distances, angles, angle_unit="deg"):
    obs = []
    for st in ["A","B","D","E"]:
        if st not in known:
            raise ValueError(f"未知測站 {st} 沒有座標")
        if st not in distances:
            raise ValueError(f"距離缺少 {st}")
        if st not in angles:
            raise ValueError(f"角度缺少 {st}")
        s = float(distances[st])
        a = float(angles[st])
        if angle_unit.lower().startswith("deg"):
            a = deg2rad(a)
        # 角度觀測就是「從參考線 i->Ref[i] 旋到 i->C」的角量測值
        obs.append({"station": st, "dist": s, "angle": a})
    return obs

# -------------------------------
# 最小二乘（等權）
# -------------------------------
def solve_unknown(known, obs, Ref, X0=None, Y0=None, max_iter=30, tol=1e-10):
    # 初值：幾何中心（若未給）
    if X0 is None or Y0 is None:
        X0 = np.mean([p[0] for p in known.values()])
        Y0 = np.mean([p[1] for p in known.values()])

    X, Y = float(X0), float(Y0)

    for it in range(max_iter):
        A_list, w_list = [], []
        for o in obs:
            st = o["station"]
            Xi, Yi = known[st]
            dX, dY = X - Xi, Y - Yi
            r2 = dX*dX + dY*dY
            r  = math.sqrt(r2)

            # --- 距離方程 ---
            # w = s_obs - r
            w_d = o["dist"] - r
            A_list.append([ dX/r,  dY/r ])
            w_list.append(w_d)

            # --- 角度方程 ---
            # 角度定義：alpha_obs = Az(i->C) - Az(i->Ref)
            # 線性化殘差用：w_th = alpha_obs - (phi - Theta)，最後 wrap 到 (-pi,pi]
            phi   = math.atan2(dY, dX)
            refX, refY = known[Ref[st]]
            Theta = azimuth(Xi, Yi, refX, refY)
            w_th  = wrap_pi(o["angle"] - (phi - Theta))
            # 偏導
            A_list.append([ dY/r2, -dX/r2 ])
            w_list.append(w_th)

        A = np.array(A_list, float)
        w = np.array(w_list,  float).reshape(-1,1)

        # 正規方程（等權）
        N = A.T @ A
        u = A.T @ w
        dx = np.linalg.solve(N, u)

        X += dx[0,0]
        Y += dx[1,0]

        if max(abs(dx[0,0]), abs(dx[1,0])) < tol:
            break

    # 殘差與精度
    v = A @ dx - w
    dof = len(w) - 2
    sigma0 = math.sqrt(float(v.T @ v) / max(dof,1))
    Cov = sigma0**2 * np.linalg.inv(N)
    sigX = math.sqrt(Cov[0,0])
    sigY = math.sqrt(Cov[1,1])
    rho  = Cov[0,1]/(sigX*sigY)

    return {"X": X, "Y": Y, "sigma0": sigma0,
            "sigma_X": sigX, "sigma_Y": sigY, "rho": rho,
            "iterations": it+1}

# -------------------------------
# 測量數據
# -------------------------------
# 假設測量距離
user_dist = {
    "A": 5.25,
    "B": 4.90,
    "D": 3.10,
    "E": 4.35,
}
# 假設角度是度（例如 A 站 ∠CAB = 38.2° ）
user_ang_deg = {
    "A": 38.2,
    "B": 121.5,
    "D": 210.7,
    "E": 296.4,
}

# 建立觀測
obs_user = build_obs_from_user(known, user_dist, user_ang_deg, angle_unit="deg")

# 給一個合理初值
X0 = 4.2
Y0 = 3.8

ans = solve_unknown(known, obs_user, Ref, X0=X0, Y0=Y0)
print(ans)
