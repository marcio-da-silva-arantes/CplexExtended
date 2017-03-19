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
import ilog.cplex.IloCplex.IncumbentCallback;

/**
 *
 * @author marcio
 */
public abstract class IncumbentBest extends IncumbentCallback implements iCplexExtract{
    public abstract void main2() throws Exception;
    
    private IncumbentBest next;
    private IncumbentBest back;
    public void setNext(IncumbentBest next) {
        this.next = next;
    }
    public void setBack(IncumbentBest back) {
        this.back = back;
    }
    private IncumbentBest orig(){
        IncumbentBest orig = this;
        while(orig.back!=null){
            orig = orig.back;
        }
        return orig;
    }
    
    @Override
    public final void main() throws IloException{
        //System.out.println("MAIN");
        try {                
            main2();
        } catch (Throwable ex) {
            ex.printStackTrace();
        }
        if(next!=null){
            next.main();
        }
    }

    @Override
    public IloCplex.Status getStatus() throws IloException {
        //return orig().getStatus();//IloCplex.Status.Feasible;
        IncumbentBest orig = orig();
        if(orig==this){
            return IloCplex.Status.Feasible;
        }else{
            return orig.getStatus(); //To change body of generated methods, choose Tools | Templates.
        }
    }
    @Override
    public double getBestObjValue() throws IloException {
        IncumbentBest orig = orig();
        if(orig==this){
            return super.getBestObjValue();
        }else{
            return orig.getBestObjValue(); //To change body of generated methods, choose Tools | Templates.
        }
    }
    @Override
    public double getObjValue() throws IloException {
        IncumbentBest orig = orig();
        if(orig==this){
            return super.getObjValue();
        }else{
            return orig.getObjValue(); //To change body of generated methods, choose Tools | Templates.
        }
    }
    
    @Override
    public double getValue(IloNumExpr var) throws IloException {
        IncumbentBest orig = orig();
        if(orig==this){
            return super.getValue(var);
        }else{
            return orig.getValue(var); //To change body of generated methods, choose Tools | Templates.
        }
    }
    @Override
    public double getValue(IloNumExpr V, boolean ignoreUnknown) throws IloException{
        if(ignoreUnknown){
            try{
                return getValue(V);
            }catch(IloCplex.UnknownObjectException ex){
                if(V instanceof IloNumVar){
                    IloNumVar v = (IloNumVar) V;
                    if(Math.abs(v.getUB()-v.getLB())<1e-4){
                        return v.getUB();
                    }
                }
                return Double.NaN;
            }catch(java.lang.NullPointerException ex){
                if(V instanceof IloNumVar){
                    IloNumVar v = (IloNumVar) V;
                    if(Math.abs(v.getUB()-v.getLB())<1e-4){
                        return v.getUB();
                    }
                }
                return Double.NaN;
            }
        }else{
            return getValue(V);
        }
    }
    @Override
    public double[] getValues(IloNumExpr[] V) throws IloException {
        return getValues(V, false);
    }
    @Override
    public double[] getValues(IloNumExpr V[], boolean ignoreUnknown) throws IloException{
        double X[] = new double[V.length];
        for(int i=0; i<V.length; i++){
            X[i] = getValue(V[i], ignoreUnknown);
        }
        return X;
    }
  
    @Override
    public double[][] getValues(IloNumExpr[][] V) throws IloException {
        return getValues(V, false);
    }
    @Override
    public double[][] getValues(IloNumExpr V[][], boolean ignoreUnknown) throws IloException{
        double X[][] = new double[V.length][];
        for(int i=0; i<V.length; i++){
            X[i] = getValues(V[i], ignoreUnknown);
        }
        return X;
    }
    
    
    @Override
    public double[][][] getValues(IloNumExpr V[][][]) throws IloException{
        return getValues(V, false);
    }
    @Override
    public double[][][] getValues(IloNumExpr V[][][], boolean ignoreUnknown) throws IloException{
        double X[][][] = new double[V.length][][];
        for(int i=0; i<V.length; i++){
            X[i] = getValues(V[i], ignoreUnknown);
        }
        return X;
    }
    @Override
    public double[][][][] getValues(IloNumExpr V[][][][]) throws IloException{
        return getValues(V, false);
    }
    @Override
    public double[][][][] getValues(IloNumExpr V[][][][], boolean ignoreUnknown) throws IloException{
        double X[][][][] = new double[V.length][][][];
        for(int i=0; i<V.length; i++){
            X[i] = getValues(V[i], ignoreUnknown);
        }
        return X;
    }
    
    @Override
    public void print(String name, String format, double[] V) {
        System.out.println("------------------------ [ "+name+"] ------------------------");
        for(double v : V){
            System.out.printf(format+" ", v);
        }
        System.out.println();
    }
    @Override
    public void print(String name, String format, double[][] V) {
        System.out.println("------------------------ [ "+name+"] ------------------------");
        for (double[] v : V) { 
            for (int j = 0; j < v.length; j++) {
                System.out.printf(format+" ", v[j]);
            }
            System.out.println();
        }
    }
   @Override
    public void print(String name, String format, double[][]... V) {
        System.out.println("------------------------ [ "+name+"] ------------------------");
        int d2 = 0;
        for (double[][] v : V) {
            d2 = Math.max(d2, v.length);
        }
        for(int t=0; t<d2; t++){
            for (int n = 0; n < V.length; n++) {
                if(t<V[n].length){
                    for (int j = 0; j < V[n][t].length; j++) {
                        System.out.printf(format+" ", V[n][t][j]);
                    }
                }else{
                    for (int j = 0; j < V[n][0].length; j++) {
                        String s = String.format(format, 0.0);
                        System.out.printf("%"+s.length()+"s ", "-");
                    }
                }
                System.out.print("| ");
            }
            System.out.println();
        }
    }
}
