# 同位素峰簇计算器

仅用 Python 标准库实现：根据离线元素同位素丰度表，计算分子式的整数总质量数概率分布（同位素峰簇），输出汇总 JSON 与逐峰 CSV。

## 用法

```
python -m isotopes 输入.json --out 输出目录
```

输出目录下生成：

- `summary.json`：分子式、电荷、合并元素计数、总概率、平均分子质量、平均 m/z 及全部峰。
- `peaks.csv`：逐峰 `mass_number,probability,mean_exact_mass,mz`，按总质量数递增。

任一非法输入都会在 stderr 报出分子或元素位置并以非零码退出，且不生成成功报告。

## 输入格式

```json
{
  "formula": "K4(ON(SO3)2)2",
  "charge": 2,
  "elements": {
    "S": [
      {"mass_number": 32, "exact_mass": 31.972071, "abundance": 0.9499},
      {"mass_number": 33, "exact_mass": 32.971459, "abundance": 0.0075}
    ]
  }
}
```

- `formula`：元素符号 + 正整数下标 + 可嵌套括号，如 `H2O`、`CH3COOH`、`K4(ON(SO3)2)2`。省略下标视为 1；分子式与括号内容不能为空，下标不得为 0；不支持小数系数与离子符号；未知元素、括号未闭合均报错。
- `charge`：正整数电荷数。
- `elements`：每个元素一个同位素数组。`mass_number` 为正整数且同元素内唯一；`exact_mass` 为有限正数（单位 u）；`abundance` 为有限非负数，每个元素丰度和与 1 的差不得超过 1e-9，计算前自动归一化。

## 合并算法

假设各原子的同位素选择相互独立。从 `{0: (概率1, 精确质量加权和0)}` 出发，对分子中每个原子与其元素的同位素分布做卷积：

- 新概率 `p' = p * abundance`；
- 新精确质量加权和 `m' = m * abundance + p' * exact_mass`；
- 以整数总质量数为键合并所有组合（相同键的概率与加权和分别累加）。

最终每个峰的概率加权精确质量均值为 `m / p`，`m/z = (m / p) / charge`。不按精确质量近似分箱，不裁剪低丰度峰，仅输出正概率峰；忽略电子质量与加合离子。汇总中的平均分子质量为所有峰的精确质量加权和（总概率归一后即为期望质量）。

## 自测

```
python -m unittest discover -s tests
```

## 示例

```
python -m isotopes examples/water.json --out out_water
```
