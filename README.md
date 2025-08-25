# 四個已知點（角度 + 距離）解未知點
![示意圖](images/diagram.png)

- **未知點**：$C(X,Y)$  
- **已知點**：$P_i(X_i,Y_i)$  
- **參考線終點**：$R_i(X_{R_i},Y_{R_i})$  
  - 在站 $i$ 量角時，角度定義為「從參考線 $i\to R_i$ 旋到視線 $i\to C$」。

**觀測量（8 個）：**
- 距離（4）：$s_A,s_B,s_D,s_E$  
- 角度（4）：$\alpha_A,\alpha_B,\alpha_D,\alpha_E$

## 觀測方程式

定義：
$$
\Delta X_i = X - X_i, \quad 
\Delta Y_i = Y - Y_i, \quad 
r_i = \sqrt{\Delta X_i^2 + \Delta Y_i^2}
$$

- 視線方位角：
$$
\phi_i = \operatorname{atan2}(\Delta Y_i,\Delta X_i)
$$

- 參考線方位角（常數）：
$$
\Theta_i = \operatorname{atan2}(Y_{R_i}-Y_i,\;X_{R_i}-X_i)
$$

**觀測方程式：**
$$
\begin{aligned}
s_i^{obs} &= r_i + v_i^{(d)} \\[6pt]
\alpha_i^{obs} &= (\phi_i - \Theta_i) + v_i^{(\theta)}
\end{aligned}
$$

---

## 3. 線性化與偏導

### 距離式
$$
\frac{\partial r_i}{\partial X} = \frac{\Delta X_i}{r_i}, \qquad
\frac{\partial r_i}{\partial Y} = \frac{\Delta Y_i}{r_i}
$$

線性化：
$$
w_i^{(d)} = s_i^{obs} - r_i(X_0,Y_0) \;\approx\; -\frac{\Delta X_i}{r_i} dX - \frac{\Delta Y_i}{r_i} dY + v_i^{(d)}
$$

### 角度式
$$
\frac{\partial \phi_i}{\partial X} = -\frac{\Delta Y_i}{r_i^2}, \qquad
\frac{\partial \phi_i}{\partial Y} = \frac{\Delta X_i}{r_i^2}
$$

線性化：
$$
w_i^{(\theta)} = \alpha_i^{obs} - (\phi_i - \Theta_i)|_{(X_0,Y_0)} \;\approx\; \frac{\Delta Y_i}{r_i^2} dX - \frac{\Delta X_i}{r_i^2} dY + v_i^{(\theta)}
$$

---

## 4. $A$ 矩陣與 $w$ 向量

$$
\mathbf{A}=
\begin{bmatrix}
-\tfrac{\Delta X_A}{r_A} & -\tfrac{\Delta Y_A}{r_A} \\
-\tfrac{\Delta X_B}{r_B} & -\tfrac{\Delta Y_B}{r_B} \\
-\tfrac{\Delta X_D}{r_D} & -\tfrac{\Delta Y_D}{r_D} \\
-\tfrac{\Delta X_E}{r_E} & -\tfrac{\Delta Y_E}{r_E} \\[6pt]
\tfrac{\Delta Y_A}{r_A^2} & -\tfrac{\Delta X_A}{r_A^2} \\
\tfrac{\Delta Y_B}{r_B^2} & -\tfrac{\Delta X_B}{r_B^2} \\
\tfrac{\Delta Y_D}{r_D^2} & -\tfrac{\Delta X_D}{r_D^2} \\
\tfrac{\Delta Y_E}{r_E^2} & -\tfrac{\Delta X_E}{r_E^2}
\end{bmatrix}
$$

$$
\mathbf{w}=
\begin{bmatrix}
s_A^{obs} - r_A \\
s_B^{obs} - r_B \\
s_D^{obs} - r_D \\
s_E^{obs} - r_E \\[6pt]
\alpha_A^{obs} - (\phi_A - \Theta_A) \\
\alpha_B^{obs} - (\phi_B - \Theta_B) \\
\alpha_D^{obs} - (\phi_D - \Theta_D) \\
\alpha_E^{obs} - (\phi_E - \Theta_E)
\end{bmatrix}
$$

---

## 5. 最小二乘解

等權情況下：
$$
\delta \mathbf{x} = (\mathbf{A}^\mathsf{T}\mathbf{A})^{-1}\mathbf{A}^\mathsf{T}\mathbf{w}, \qquad
\mathbf{x} \leftarrow \mathbf{x} + \delta\mathbf{x}
$$

---

## 6. 精度評估

- 自由度：$f = 8 - 2 = 6$  
- 單位權中誤差：
$$
\hat\sigma_0 = \sqrt{\frac{\mathbf{e}^\mathsf{T}\mathbf{e}}{f}}
$$
- 協方差矩陣：
$$
\Sigma_{\hat x} = \hat\sigma_0^2 (\mathbf{A}^\mathsf{T}\mathbf{A})^{-1}
$$

---