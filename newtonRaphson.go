package main

import (
    "errors"
    "fmt"
    "math"
)

//牛顿拉弗森迭代算法解非线性方程
// xn = x - f(x)/f'(x) 迭代方程

func calcY(factors []float64, x float64, order int) (ret float64) {
    ret = 0
    for i := 0; i <= order; i++ {
        ret += factors[i] * math.Pow(x, float64(order - i))
    }
    return ret
}
// newtonRaphson 牛顿哈夫森算法解一元多次方程
func newtonRaphson(factors []float64, order int, x0 float64, precision float64) (float64, error) {
    diffFactors := make([]float64, order)
    for i := 0; i < order; i++ {
        diffFactors[i] = factors[i] * float64(order - i)
    }
    y0 := calcY(factors, x0, order)
    diffy0 := calcY(diffFactors, x0, order - 1)
    var xn, yn, temp, diffyn float64
    xn = x0 - y0 / diffy0
    temp = x0
    var counter = 1
    for {
        if math.Abs(xn - temp) < precision {
            break
        }
        temp = xn
        yn = calcY(factors, xn, order)
        diffyn = calcY(diffFactors, xn, order - 1)
        xn = xn - yn / diffyn
        counter++
        if counter > 1000 {
            return math.NaN(), errors.New("无法收敛，此点附近无解 ")
        }
    }

    return xn, nil
}

func main() {
    // example
    factors  := []float64{1, 1, -2, 1, 2, 0}
    x0 := -1.65
    precision := 1e-15
    ret, err := newtonRaphson(factors, 5, x0, precision)
    if err != nil {
        fmt.Println(err)
    }
    fmt.Println("方程解为：", ret) //

    // 验证方程解
    yret := calcY(factors, ret,5) // 预期0
    fmt.Println("方程解对应的y: ", yret)
}
