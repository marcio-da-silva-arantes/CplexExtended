/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package CplexExtended;

import ilog.concert.IloException;
import ilog.concert.IloNumExpr;
import ilog.concert.IloNumVar;


/**
 *
 * @author marcio
 */
public class Hyperplane {
    public final double a[];
    public final double b;
    public Hyperplane(double b, double ...a) {
        this.a = a;
        this.b = b;
    }
    public double scalarProd(double ...p) {
        double sum = 0;
        int length = Math.min(a.length, p.length);
        for(int n=0; n<length; n++){
            sum += a[n] * p[n];
        }
        return sum;
    }

    public IloNumExpr scalProd(CplexExtended cplex, IloNumVar... x) throws IloException {
        IloNumExpr sum = null;
        int length = Math.min(a.length, x.length);
        for(int n=0; n<length; n++){
            sum = cplex.SumProd(sum, a[n], x[n]);
        }
        return sum;
    }
    
}
