/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package CplexExtended;

import ilog.concert.IloException;
import ilog.concert.IloNumExpr;
import ilog.concert.IloNumVar;
import ilog.cplex.IloCplex;
import ilog.cplex.IloCplex.Status;

/**
 *
 * @author marcio
 */
public interface iCplexExtract {
    public double getBestObjValue() throws IloException;
    public double getObjValue() throws IloException;
    public double getValue(IloNumExpr expr) throws IloException;
    //public double getValue(IloNumVar var) throws IloException;
    public double getValue(IloNumExpr var, boolean ignoreUnknown) throws IloException;
    public double[] getValues(IloNumExpr V[]) throws IloException;
    public double[] getValues(IloNumExpr V[], boolean ignoreUnknown) throws IloException;
    public double[][] getValues(IloNumExpr V[][]) throws IloException;
    public double[][] getValues(IloNumExpr V[][], boolean ignoreUnknown) throws IloException;
    public double[][][] getValues(IloNumExpr V[][][]) throws IloException;
    public double[][][] getValues(IloNumExpr V[][][], boolean ignoreUnknown) throws IloException;
    public double[][][][] getValues(IloNumExpr V[][][][]) throws IloException;
    public double[][][][] getValues(IloNumExpr V[][][][], boolean ignoreUnknown) throws IloException;
    
    public void print(String name, String format, double[] V);
    public void print(String name, String format, double[][] V);
    public void print(String name, String format, double[][] ...V);
    public Status getStatus() throws IloException;
}
