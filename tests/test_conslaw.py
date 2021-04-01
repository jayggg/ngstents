import sys

sys.path.append("../demo/wave")
sys.path.append("../demo/advection")
sys.path.append("../demo/symbolic")

def test_wave2d():
    ''' 
    reference values generated with NGSolve-6.2.2102-17-g02efd4be3
    and ngstents commit 55862b9
    
    standing wave example:
    * mesh size h = 0.25, spatial order = 2
    * structure-aware Taylor time stepping with 3 stages and 8 substeps within each tent
    * reference for spatial L2 error: 3.60383e-04
    '''
    import wave2d
    assert wave2d.l2error <= 4e-4

def test_wave2d_timdepbc():
    '''
    plane wave example:
    * mesh size h = 0.25, spatial order = 2
    * structure-aware Taylor time stepping with 3 stages and 8 substeps within each tent
    * reference for spatial L2 error: 1.11393e-02
    '''
    import wave2d_timedepbc
    assert wave2d_timedepbc.l2error <= 2e-2
    

def test_advection2d():
    ''' 
    reference values generated with NGSolve-6.2.2102-17-g02efd4be3
    and ngstents commit 55862b9
    
    standing wave example:
    * mesh size h = 0.15, spatial order = 4
    * structure-aware Taylor time stepping with 5 stages and 8 substeps within each tent
    * reference for spatial L2 error: 2.79097e-03
    '''
    import advection2d
    assert advection2d.l2error <= 3e-3

def test_symbolic_wave():
    ''' 
    reference values generated with NGSolve-6.2.2102-17-g02efd4be3
    and ngstents commit 55862b9
    
    standing wave example:
    * mesh size h = 0.25, spatial order = 2
    * structure-aware Runge-Kutta time stepping with 8 substeps within each tent
    * reference for spatial L2 error: 3.78149e-04
    '''
    import symbolic_wave
    assert symbolic_wave.l2error <= 4e-4

if __name__ == "__main__":
    functions = [test_wave2d, test_wave2d_timdepbc,
                 test_advection2d, test_symbolic_wave]
    passed = []
    print("Test wave equation:")
    for func in functions:
        try:
            func()
            passed.append(True)
        except AssertionError:
            passed.append(False)
    print("Function name \t\t pass/fail")
    for i in range(len(functions)):
        if passed[i]:
            print(functions[i].__name__,": \t\t pass")
        else:
            print(functions[i].__name__,": \t\t fail")
