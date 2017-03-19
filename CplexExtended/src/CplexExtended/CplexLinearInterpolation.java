/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package CplexExtended;

import ilog.concert.IloException;
import ilog.concert.IloIntVar;
import ilog.concert.IloNumExpr;
import ilog.concert.IloNumVar;

/**
 *
 * @author marci
 */
public abstract class CplexLinearInterpolation implements iCplexFunction1v{
    private final IloIntVar y[];    //model points the fuction
    private final IloNumVar w[];    //model linear interpolation bettewen points of the function
    public final IloNumExpr f;      //value returned for the function
    public final int N;
    public final String name;
    
    public CplexLinearInterpolation(CplexExtended cplex, IloNumExpr x, int N, String name) throws IloException {
        this.N = N;
        this.name = name;
        y = cplex.boolVarArray(N, name+".y");
        w = cplex.numVarArray(N, 0.0, 1.0, name+".y");
        
        // sum{y[n]} = 1
        IloNumExpr one = null;
        for(int n=0; n<N; n++){
            one = cplex.SumProd(one, 1.0, y[n]);
        }
        cplex.addEq(one, 1, name+".sum(y)_1");
        
        // w[n] <= y[n] 
        for(int n=0; n<N; n++){
            cplex.addLe(w[n], y[n], name+".w(n)-y(n)");
        }
        
        // x = sum{ x(n) * y[n] } + sum { (x(n+1) - x(n)) * w[n] }
        IloNumExpr var = null;
        for(int n=0; n<N; n++){
            var = cplex.SumProd(var, x(n), y[n]);
        }
        for(int n=0; n<N; n++){
            var = cplex.SumProd(var, x(n+1)-x(n), w[n]);
        }
        cplex.addEq(x, var, name+".var");
        
        // f = sum{ f(n) * y[n] } + sum { (f(n+1) - f(n)) * w[n] }
        IloNumExpr func = null;
        for(int n=0; n<N; n++){
            func = cplex.SumProd(func, f(n), y[n]);
        }
        for(int n=0; n<N; n++){
            func = cplex.SumProd(func, f(n+1)-f(n), w[n]);
        }
        f = func;
    }
}
