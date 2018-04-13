from eddlib.estimate import golden_section_search
import numpy as np

def f1(x):
    return -1*(x - 3)**2 + 100

def f2(x):
    return -1*(x - 20)**2 + 100

def f3(x):
    return -1*(x - 35)**2 + 100

xs = np.linspace(0, 40, 10000)

def test_f1():
    idx = f1(xs).argmax()    
    correct_val = xs[idx]
    x = golden_section_search(f1, 0, 30, 40, 0.01)
    assert abs(correct_val - x) < 0.01

def test_f2():
    idx = f2(xs).argmax()
    correct_val = xs[idx]
    x = golden_section_search(f2, 0, 30, 40, 0.005)
    msg = 'correct %.2f, but got %.2f' % (correct_val, x)
    assert abs(correct_val - x) < 0.01, msg

def test_f3():
    idx = f3(xs).argmax()    
    correct_val = xs[idx]
    x = golden_section_search(f3, 0, 30, 40, 0.01)
    assert abs(correct_val - x) < 0.01
    
