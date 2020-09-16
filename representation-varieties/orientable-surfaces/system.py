import sympy as sp

q = sp.Symbol('q')

# Util functions
def eq_subs_by_fraction(eq, x, a, b):
    deg = sp.degree(eq, x)
    return sp.simplify(eq.subs(x, a / b) * (b ** deg))

def is_square(t):
    t = sp.factor(t)
    if t.is_complex:
        return sp.sqrt(t)
    if t.is_Pow and t.args[1] % 2 == 0:
        return sp.Pow(t.args[0], t.args[1] / 2)
    if t.is_Mul:
        sqrt_t = 1
        for f in t.args:
            sqrt_f = is_square(f)
            if sqrt_f == None:
                return None
            sqrt_t *= sqrt_f
        return sqrt_t
    return None

# System class
class System:
    # variables is a set of sympy symbols used in the equations
    # cl_eqs is a list of sympy expressions that should equal zero
    # op_eqs is a list of sympy expressions that should NOT equal zero
    
    def __init__(self, variables = {}, cl_eqs = [], op_eqs = []):
        self.vars = set(variables)
        self.cl_eqs = list(cl_eqs) 
        self.op_eqs = list(op_eqs)
    
    def check_unsatisfiable(self):
        # Check for equations of the form 'non-zero number = 0'
        for eq in [ e for e in self.cl_eqs if e.is_complex and e != 0 ]:
            return 0
        
        # Check for equations of the form '0 != 0'
        if 0 in self.op_eqs:
            return 0
                
        return None
    
    def check_unused_variables(self): # TODO: change to: check if system of equation can be factored into S1 x S2 x ... x Sn?
        # If there are variables that do not appear in any equation, factor them out
        used_symbols = set.union(set(), *[ eq.free_symbols for eq in self.cl_eqs + self.op_eqs ])
        unused_symbols = self.vars.difference(used_symbols)
        if len(unused_symbols) > 0:
            self.vars = self.vars.difference(unused_symbols)
            return sp.expand((q ** len(unused_symbols)) * self.compute_class())

        return None
        
    def check_filter_equations(self):
        made_changes = False
        
        # Expand, and remove instances of '0 = 0' and 'non-zero number != 0'
        len_cl, len_op = len(self.cl_eqs), len(self.op_eqs)
        self.cl_eqs = list({ sp.expand(e) for e in self.cl_eqs if e != 0 })
        self.op_eqs = list({ sp.expand(e) for e in self.op_eqs if not e.is_complex or e == 0 })
        if len(self.cl_eqs) != len_cl or len(self.op_eqs) != len_op:
            made_changes = True
        
        # If any equation is of the form 'x = 0', we can solve for x = 0 and continue with the remaining system
        for x in [ eq for eq in self.cl_eqs if eq.is_Symbol and eq in self.vars ]:
            self.vars = self.vars.difference({x,})
            self.cl_eqs = [ e.subs(x, 0) for e in self.cl_eqs ]
            self.op_eqs = [ e.subs(x, 0) for e in self.op_eqs ]
            made_changes = True
        
        for eq in self.cl_eqs:
            eq_factor = sp.factor(eq)
        
            # If any equation is of the form '(...)^n = 0', we can simply substitute with '(...) = 0' (provided n is positive)
            if eq_factor.is_Pow:
                if eq_factor.args[1] < 0:
                    return 0
                
                self.cl_eqs = [ eq_factor.args[0] if e == eq else e for e in self.cl_eqs]
            
            # If any equation is a product, we can remove all denominators and all numeric factors. Also, one can remove all higher multiplicities
            if eq_factor.is_Mul:
                unnecessary_factors = [ a for a in eq_factor.args if a.is_complex or (a.is_Pow and a.args[1] < 0) ]
                higher_multiplicity = [ a for a in eq_factor.args if a.is_Pow and a.args[1] > 0 ]
                if unnecessary_factors or higher_multiplicity:
                    replacement = sp.Mul(*( f.args[0] if f in higher_multiplicity else f for f in eq_factor.args if f not in unnecessary_factors ))
                    self.cl_eqs = [ replacement if e == eq else e for e in self.cl_eqs ]
                    made_changes = True
        
        for eq in self.op_eqs:
            eq_factor = sp.factor(eq)
            
            # If an equation is of the form '(u) * ... * (v) != 0', replace by '(u) != 0 & ... & (v) != 0'
            if eq_factor.is_Mul:
                self.op_eqs = [ e for e in self.op_eqs if e != eq ] + [ f for f in eq_factor.args if not f.is_complex]
                made_changes = True
        
        if made_changes:
            return self.compute_class()
        
        return None
    
    def check_if_point(self):
        # Check if system describes a point (i.e. no variables and no equations)
        if not self.cl_eqs and not self.op_eqs and not self.vars:
            return 1

        return None
    
    def check_single_variable_equations(self):
        # Check for (closed) equations with only one free variable, then solve for x
        for eq in self.cl_eqs:
            symbols = eq.free_symbols
            if len(symbols) != 1:
                continue
            
            x = symbols.pop()
            if x not in self.vars:
                continue

            # Simply try all solutions and add the classes
            solutions_for_x = sp.solve(eq, x)
            total = 0
            succeeded = True
            for v in solutions_for_x:
                c = System(self.vars.difference({x,}), [ e.subs(x, v) for e in self.cl_eqs if e != eq ], [ e.subs(x, v) for e in self.op_eqs ]).compute_class()
                if c == None: # if any fails, break
                    succeeded = False
                    break
                total += c
            if succeeded:
                return total
            
        return None
            
    def check_linear_equations(self):
        # Look for something of the form 'x * (w) + (u) = 0' with x not in (w), (u)
        for eq in self.cl_eqs:
            eq_expand = sp.expand(eq)
            if not eq_expand.is_Add:
                continue

            for x in [ s for s in eq_expand.free_symbols if s in self.vars and sp.degree(eq_expand, s) == 1]:
                u = sp.Add(*[ t for t in eq_expand.args if x not in t.free_symbols ])
                w = sp.expand((eq_expand - u) / x)
                
                # Now either (w) = 0 (implying (u) = 0 as well):
                if w.is_complex:
                    c_1 = 0
                else:
                    c_1 = System(self.vars, [ e for e in self.cl_eqs if e != eq ] + [ w, u ], self.op_eqs).compute_class()
                    if c_1 == None:
                        continue

                # Or (w) != 0 and we can use the equation to solve for x = -u / w:
                c_2_cl_eqs = [ eq_subs_by_fraction(e, x, -u, w) for e in self.cl_eqs if e != eq ]
                c_2_op_eqs = [ eq_subs_by_fraction(e, x, -u, w) for e in self.op_eqs ] + [ w ]
                c_2 = System(self.vars.difference({x,}), c_2_cl_eqs, c_2_op_eqs).compute_class()
                if c_2 == None:
                    continue
                return c_1 + c_2
            
        return None
    
    def check_product_equations(self):
        # Check for equations of the form '(v) * (w1) * ... * (wn) = 0'
        for eq in self.cl_eqs:
#             eq_factor = sp.factor(eq, extension=[sp.I])
            eq_factor = sp.factor(eq)
            if not eq_factor.is_Mul:
                continue
            
            v = eq_factor.args[0]
            prod_wi = sp.cancel(eq_factor / v)
                        
            # Either (v) = 0:
            c_1 = System(self.vars, [ e for e in self.cl_eqs if e != eq ] + [ v ], self.op_eqs).compute_class()
            if c_1 == None:
                continue
            # Or (v) != 0 and (prod_wi) = 0
            c_2 = System(self.vars, [ e for e in self.cl_eqs if e != eq] + [ prod_wi ], self.op_eqs + [ v ]).compute_class()
            if c_2 == None:
                continue
            return c_1 + c_2
        
        return None
    
    def check_quadratic_equations(self):
        # Look for something quadratic to solve 'x^2 * (u) + x * (v) + (w) = 0'
        for eq in self.cl_eqs:
            eq_expand = sp.expand(eq)
            if not eq_expand.is_Add:
                continue

            for x in [ s for s in eq_expand.free_symbols if s in self.vars and sp.degree(eq_expand, s) == 2]:                
                w = sp.Add(*[ t for t in eq_expand.args if x not in t.free_symbols ])
                v = sp.cancel(sp.Add(*[ t for t in eq_expand.args if sp.degree(t, x) == 1 ]) / x)
                u = sp.cancel((eq_expand - w - v * x) / x**2)
                
                # The discriminant must be a square, otherwise the isomorphism we want to setup is not an algebraic map!
                discr = v**2 - 4*u*w
                sqrt_D = is_square(discr)
                if sqrt_D == None:
                    continue
                
                # Case (u) = 0:
                if u.is_complex:
                    c_u_is_zero = 0
                else:
                    c_u_is_zero = System(self.vars, [ e for e in self.cl_eqs if e != eq ] + [ u, x * v + w ], self.op_eqs).compute_class()
                    if c_u_is_zero == None:
                        continue
                
                # Case (u) != 0:
                # Case A: discr = 0 and x = -v / (2 * u)
                c_A_cl_eqs = [ eq_subs_by_fraction(e, x, -v, 2*u) for e in self.cl_eqs if e != eq ] + [ discr ]
                c_A_op_eqs = [ eq_subs_by_fraction(e, x, -v, 2*u) for e in self.op_eqs ] + [ u ]
                c_A = System(self.vars.difference({x,}), c_A_cl_eqs, c_A_op_eqs).compute_class()
                if c_A == None:
                    continue
                # Case B: discr != 0 and x = (-v + sqrt_D) / (2 * u)
                c_B_cl_eqs = [ eq_subs_by_fraction(e, x, -v + sqrt_D, 2*u) for e in self.cl_eqs if e != eq ]
                c_B_op_eqs = [ eq_subs_by_fraction(e, x, -v + sqrt_D, 2*u) for e in self.op_eqs ] + [ u, discr ]
                c_B = System(self.vars.difference({x,}), c_B_cl_eqs, c_B_op_eqs).compute_class()
                if c_B == None:
                    continue
                # Case C: discr != 0 and x = (-v - sqrt_D) / (2 * u)
                c_C_cl_eqs = [ eq_subs_by_fraction(e, x, -v - sqrt_D, 2*u) for e in self.cl_eqs if e != eq ]
                c_C_op_eqs = [ eq_subs_by_fraction(e, x, -v - sqrt_D, 2*u) for e in self.op_eqs ] + [ u, discr ]
                c_C = System(self.vars.difference({x,}), c_C_cl_eqs, c_C_op_eqs).compute_class()
                if c_C == None:
                    continue
                    
                return c_u_is_zero + c_A + c_B + c_C
            
        return None
    
    def check_open_equations(self):
        #                              case B          case A
        # Using that { x != 0 } = { x is whatever } - { x = 0}:
        for eq in self.op_eqs:
            # Case A: eq = 0
            c_A = System(self.vars, self.cl_eqs + [ eq ], [ e for e in self.op_eqs if e != eq ]).compute_class()
            if c_A == None:
                continue
            # Case B: eq = whatever
            c_B = System(self.vars, self.cl_eqs, [ e for e in self.op_eqs if e != eq ]).compute_class()
            if c_B == None:
                continue
                
            return c_B - c_A
        
        return None
    
    def compute_class(self):        
#         hash_ = '{' + (','.join([str(e) + ' = 0' for e in self.cl_eqs] + [str(e) + ' != 0' for e in self.op_eqs])) + '} in ' + str(self.vars) + ''
#         print(hash_)
        
        c = self.check_unsatisfiable()
        if c != None:
            return c
    
        c = self.check_if_point()
        if c != None:
            return c

        c = self.check_unused_variables()
        if c != None:
            return c
        
        c = self.check_filter_equations()
        if c != None:
            return c
                
        c = self.check_single_variable_equations()
        if c != None:
            return c

        c = self.check_linear_equations()
        if c != None:
            return c
        
        c = self.check_product_equations()
        if c != None:
            return c

        c = self.check_open_equations()
        if c != None:
            return c
        
        return make_symbol_from_unsolved(self)
#         return None

symbols_for_unsolved = {}

def make_symbol_from_unsolved(sys):
    hash_ = '{' + (','.join([str(e) + ' = 0' for e in sys.cl_eqs] + [str(e) + ' != 0' for e in sys.op_eqs])) + '} in ' + str(sys.vars) + ''
    
    if hash_ in symbols_for_unsolved:
        return symbols_for_unsolved[hash_]
    else:    
        i = len(symbols_for_unsolved)
        x = sp.Symbol('x_' + str(i))
        symbols_for_unsolved[hash_] = x
        return x

    
    
# Compute cases
def compute_cases(variables, cl_eqs, op_eqs, cases):
    n = len(cases)
    total = 0
    for i in range(n):
        case = cases[i]
        extra_cl = case[0]
        extra_op = case[1]
        weight = case[2]
        print('Case {}/{}: ... '.format(i + 1, n), end = '')
        c = System(variables, cl_eqs + extra_cl, op_eqs + extra_op).compute_class()
        print(str(sp.factor(c)))
        total += weight * c
    return total