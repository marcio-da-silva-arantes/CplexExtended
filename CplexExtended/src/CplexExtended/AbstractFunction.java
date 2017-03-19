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
 * @author marci
 */
public abstract class AbstractFunction {
    public final String name;
    public final double lb;
    public final double ub;
    
    /**
     * A abstract representation for a function: f(x)
     * @param name  a string name for the function
     * @param lb    the minimum value for the domain 
     * @param ub    the maximum value for the domain
     */
    public AbstractFunction(String name, double lb, double ub) {
        this.name = name;
        this.lb = lb;
        this.ub = ub;
    }
    public abstract double f(double x);
    
    
    public static void main(String[] args) throws IloException {
        //testFunction(args);
        
        testComplexity(args);
        
    }
    public static void testFunction(String[] args) throws IloException {
        
        
        //defining a function: f(x) = x*sin(10*PI*x) + 1 on domain x in [-1 , +2]
        AbstractFunction func = new AbstractFunction("func", -1.0, +2.0) {
            @Override
            public double f(double x) {
                return x*Math.sin(10*Math.PI*x)+1;
            }
        };
        
        //using cplex to solve the maximum value for this the function
        CplexExtended cplex = new CplexExtended();
       
        //delcare the variable on domain lb = -1.0 and ub = +2.0
        IloNumVar x = cplex.numVar(-1, 2, "x");
        
        //create the expression for this function using 100000 linear segments
        IloNumExpr obj = cplex.Function(x, 12500, func);
        
        //added to maximize this expression
        cplex.addMaximize(obj);
        
        
        if(cplex.solve()){
            cplex.setOut(null);
            cplex.setWarning(null);
            System.out.println("Solved, status  = "+cplex.getStatus());
            System.out.println("Objective       = "+cplex.getObjValue());
            System.out.println("x*              = "+cplex.getValue(x));
            System.out.println("f*              = "+cplex.getValue(obj));
            
//            double error = 0.0015;
//            double r = 0.000;
//            for(int n=0; n<func.N; n++){
//                x.setLB(func.x(n)+r);
//                x.setUB(func.x(n)+r);
//                cplex.solve();
//                System.out.printf("sin(%g) = %g | model = %g\n", func.x(n)+r, FUNC(func.x(n)+r), cplex.getValue(func.f));
//                error = Math.max(error, Math.abs(FUNC(func.x(n)+r) - cplex.getValue(func.f)));
//            }
//            System.out.println("Maximun error = "+error);
        }else{
            System.out.println("Not solved, status = "+cplex.getStatus());
        }
    }
    
    public static void testComplexity(String[] args) throws IloException {
        //defining a function: f(x) = x*sin(10*PI*x) + 1 on domain x in [-1 , +2]
        AbstractFunction func = new AbstractFunction("func", -1.0, +2.0) {
            @Override
            public double f(double x) {
                return x*Math.sin(10*Math.PI*x)+1;
            }
        };
        
        //using cplex to solve the maximum value for this the function
        CplexExtended cplex = new CplexExtended();
        cplex.setOut(null);
        cplex.setOut(null);
        
        System.out.printf("%10s %10s %10s %10s %10s\n", "N", "Time", "X*", "F*", "factor");
        double time_old = 0;
        for(int n=0; n<8; n++){
            int N = (int)(1000*Math.pow(2, n)+0.5);
            cplex.clearModel();
            //delcare the variable on domain lb = -1.0 and ub = +2.0
            IloNumVar x = cplex.numVar(-1, 2, "x");

            //create the expression for this function using 100000 linear segments
            IloNumExpr obj = cplex.Function(x, N, func);

            //added to maximize this expression
            cplex.addMaximize(obj);
            
            
            double time = System.currentTimeMillis()/1000.0;
            if(cplex.solve()){
                time = System.currentTimeMillis()/1000.0 - time;
                System.out.printf("%10d %10g %10g %10g %10g\n", N, time, cplex.getValue(x), cplex.getObjValue(), time/time_old);
            }else{
                System.out.println("Not solved, status = "+cplex.getStatus());
            }
            time_old = time;
        }
    }
}
