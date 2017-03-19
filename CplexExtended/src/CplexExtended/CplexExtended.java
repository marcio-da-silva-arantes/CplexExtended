/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package CplexExtended;

import ilog.concert.IloException;
import ilog.concert.IloIntVar;
import ilog.concert.IloNumExpr;
import ilog.concert.IloNumVar;
import ilog.cplex.IloCplex;
import java.util.LinkedList;
import jsc.distributions.Normal;

/**
 *
 * @author marcio
 */
public class CplexExtended extends IloCplex implements iCplexExtract{
    private final double bigM;// = 1e5;
    private final double epsilon;// = 1e-4;
    
    public CplexExtended(double bigM, double epsilon) throws IloException {
        super();
        this.bigM = bigM;
        this.epsilon = epsilon;
    }
    public CplexExtended() throws IloException {
        this(1e+5, 1e-3);
    }

    
    private final LinkedList<Class> list_class = new LinkedList<Class>();
    private final LinkedList<Callback> list_cb = new LinkedList<Callback>();
    @Override
    public final void use(Callback cb) throws IloException {
        Class c = cb.getClass();
        while(c!=null && !c.getPackage().getName().equals("ilog.cplex")){
            //System.out.println(c);
            c = c.getSuperclass();
        }
        //System.out.println(c);
        if(list_class.contains(c)){
            if(cb instanceof IncumbentCallback){
                IncumbentBest old = (IncumbentBest) list_cb.get(list_class.indexOf(c));
                IncumbentBest now = (IncumbentBest) cb;
                old.setNext(now);
                now.setBack(old);
//                now.setNext(old);
//                list_cb.set(list_class.indexOf(c), now);
            }else{
                throw new IloException("Callbac type '"+c+"' override");
            }
        }else{
            super.use(cb); //To change body of generated methods, choose Tools | Templates.
            list_class.addLast(c);
            list_cb.addLast(cb);
            //System.out.println("LIST = "+list_class+", callback = "+c);
        }
    }
    
    public final void addMinimizeExtra(IloNumExpr exp) throws IloException {
        getObjective().setExpr(sum(exp,getObjective().getExpr()));
    }
    public final IloNumExpr RiskAllocation(IloNumVar delta, final double uncertainty, final double max, final int N, String name) throws IloException{
        return RiskAllocation(delta, uncertainty, max, 2.0, N, name);
    }
    /**
     * Risk allocation 
     * c(delta) = erf-inv(1-2delta) * uncertainty
     */
    public final IloNumExpr RiskAllocation(IloNumVar delta, final double uncertainty, final double max, final double base, final int N, String name) throws IloException{
        if(uncertainty<1e-6){
            return constant(0.0);
        }else{
            IloNumExpr erf_inv = ERF_Inv(delta, max, base, N, name);
            return prod(erf_inv, uncertainty);
        }
    }
    /**
     * Risk allocation 
     * c(delta) = erf-inv(1-2delta) * uncertainty
     */
    public final IloNumExpr RiskAllocation(IloNumExpr erf_inv, final double uncertainty) throws IloException{
        if(uncertainty<1e-6){
            return constant(0.0);
        }else{
            return prod(erf_inv, uncertainty);
        }
    }
    
    /**
     * erf-inv(1-2delta)
     */
    public final IloNumExpr ERF_Inv(IloNumVar delta, final double max, final double base, final int N, String name) throws IloException{
        return ERF_Inv(delta, max, base, N, true, name);
    }
    /**
     * erf-inv(1-2delta)
     */
    public final IloNumExpr ERF_Inv(IloNumVar delta, final double max, final double base, final int N, boolean bounded, String name) throws IloException{
        CplexNumFunction erf_inv = new CplexNumFunction(this, N, bounded, name) {
            @Override
            public double x(int n) {
                //return 1.0 - max/Math.pow(base, n);
                return 1.0 - 2*max/Math.pow(base, n);
            }
            @Override
            public double f(int n) {
                //System.out.println("F(x("+n+")): x = "+x(0));
                return Normal.inverseStandardCdf((1+x(n))/2.0)/Math.sqrt(2);
            }
        };
        //addLe(1, sum(prod(2,delta), erf_inv.x), name+".x_2delta");
        addEq(1, sum(prod(2,delta), erf_inv.x), name+".x_2delta");
        return erf_inv.f;
    }
    /**
     * It's used N binary variables to encode this fuction, but with a linear complexity to solve it.
     * @param alpha the angle variable in radians
     * @param lb    the minimal value for alpha
     * @param ub    the maximun value for alpha 
     * @param N     the number of segments to aproximate the function 
     * @return      a aproximation for the math function sin(alpha) 
     * @throws ilog.concert.IloException 
     */
    public final IloNumExpr sin(IloNumExpr alpha, double lb, double ub, int N) throws IloException{
        AbstractFunction sin = new AbstractFunction("sin", lb, ub) {
            @Override
            public double f(double x) {
                return Math.sin(x);
            }
        };
        return Function(alpha, N, sin);
    }
    /**
     * It's used N binary variables to encode this fuction, but with a linear complexity to solve it.
     * @param alpha the angle variable in radians
     * @param lb    the minimal value for alpha
     * @param ub    the maximun value for alpha 
     * @param N     the number of segments to aproximate the function 
     * @return      a aproximation for the math function sin(alpha) 
     * @throws ilog.concert.IloException 
     */
    public final IloNumExpr cos(IloNumExpr alpha, double lb, double ub, int N) throws IloException{
        AbstractFunction cos = new AbstractFunction("cos", lb, ub) {
            @Override
            public double f(double x) {
                return Math.cos(x);
            }
        };
        return Function(alpha, N, cos);
    }
    /**
     * It's used N binary variables to encode this fuction, but with a linear complexity to solve it.
     * @param x     the variable, domain of the function
     * @param N     the number of segments to aproximate the function 
     * @param func  a abstrat interface to represent this function
     * @return      a cplex expression that aproximate function values 
     * @throws ilog.concert.IloException 
     */
    public final IloNumExpr Function(IloNumExpr x, int N, AbstractFunction func) throws IloException{
        CplexLinearInterpolation inter = new CplexLinearInterpolation(this, x, N, func.name){
            @Override
            public double x(int n) {
                return func.lb + n*(func.ub-func.lb)/N;
            }
            @Override
            public double f(int n) {
                return func.f(x(n));
            }
        };
        return inter.f;
    }
    
    public final IloNumExpr RiskAllocationCond(IloNumVar delta, final double uncertainty, final double max, final int N, IloNumVar y, String name) throws IloException{
        return RiskAllocationCond(delta, uncertainty, max, 2.0, N, y, name);
    }
    public final IloNumExpr RiskAllocationCond(IloNumVar delta, final double uncertainty, final double max, final double base, final int N, IloNumVar y, String name) throws IloException{
        if(uncertainty<1e-6){
            return constant(0.0);
        }else{
            CplexNumFunctionCond erf_inv = new CplexNumFunctionCond(this, N, y, true, name) {
                @Override
                public double x(int n) {
                    //return 1.0 - max/Math.pow(base, n);
                    return 1.0 - 2*max/Math.pow(base, n);
                }
                @Override
                public double f(int n) {
                    //System.out.println("F(x("+n+")): x = "+x(0));
                    return Normal.inverseStandardCdf((1+x(n))/2.0)/Math.sqrt(2);
                }
            };
            //addLe(1, sum(prod(2,delta), erf_inv.x), name+".x_2delta");
            addEq(1, sum(prod(2,delta), erf_inv.x), name+".x_2delta");

            return prod(erf_inv.f, uncertainty);
        }
    }
    
    public IloNumExpr SumProd(IloNumExpr sum, double coef, IloNumExpr exp) throws IloException{
        if(exp==null){
            return sum;
        }else if(sum==null){
            return prod(coef, exp);
        }else{
            return sum(sum, prod(coef, exp));
        }
    }
    public IloNumExpr SumProd(double coef, IloNumExpr ...exp) throws IloException{
        return SumProd(null, coef, exp);
    }
    public IloNumExpr SumProd(IloNumExpr sum, double coef, IloNumExpr ...exp) throws IloException{
        for (IloNumExpr e : exp) {
            sum = SumProd(sum, coef, e);
        }
        return sum;
    }
    public IloNumExpr SumNumScalProd(IloNumExpr sum, String name, int N, double MAX, IloNumExpr... exp) throws IloException {
        return SumNumScalProd(sum, name, N, MAX, 1.0, exp);
    }
    public IloNumExpr SumNumScalProd(IloNumExpr sum, String name, int N, double MAX, double coef, IloNumExpr... exp) throws IloException {
        CplexNumScalProd func = new CplexNumScalProd(this, name, N, MAX, coef, exp);
        if(sum==null){
            return func.f;
        }else{
            return sum(sum, func.f);
        }
    }
    public IloNumExpr SumNumNorm2_2D(IloNumExpr sum, String name, int N, double MAX, IloNumExpr expX, IloNumExpr expY) throws IloException {
        return SumNumNorm2_2D(sum, name, N, MAX, 1.0, expX, expY);
    }
    public IloNumExpr SumNumNorm2_2D(IloNumExpr sum, String name, int N, double MAX, double coef, IloNumExpr expX, IloNumExpr expY) throws IloException {
        CplexNumNorm2_2D func = new CplexNumNorm2_2D(this, name, N, MAX, coef, expX, expY);
        if(sum==null){
            return func.f;
        }else{
            return sum(sum, func.f);
        }
    }
    public IloNumExpr SumNumNorm1(IloNumExpr sum, String name, double MAX, IloNumExpr... exp) throws IloException {
        return SumNumNorm1(sum, name, MAX, 1.0, exp);
    }
    public IloNumExpr SumNumNorm1(IloNumExpr sum, String name, double MAX, double coef, IloNumExpr... exp) throws IloException {
        CplexNumNorm1 func = new CplexNumNorm1(this, name, MAX, coef, exp);
        if(sum==null){
            return func.f;
        }else{
            return sum(sum, func.f);
        }
    }
    public IloNumExpr SumNumFunction(IloNumExpr sum, String name, int N, final iCplexFunction1v func, IloNumExpr... exp) throws IloException {
        for(IloNumExpr e : exp){
            CplexNumFunction temp = new CplexNumFunction(this, N, false, name.replace("(.)", "("+e+")")){
                @Override
                public double x(int n) {
                    return func.x(n);
                }
                @Override
                public double f(int n) {
                    return func.f(n);
                }
            };
            addEq(temp.x, e, name+".x_e");
            if(sum==null){
                sum = temp.f;
            }else{
                sum = sum(sum, temp.f);
            }
        }
        return sum;  
    }
    
    
    
    public static String Index(int i, int size){
        String s = "";
        int len = 0;
        while(size>0){
            size = size/10;
            len++;
        }
        size = i+1;
        while(size>0){
            size = size/10;
            len--;
        }
        while(len>0){
            s += "0";
            len--;
        }
        return s+(i+1);
    }
    /**
     * Define a array with free variables: bounds [-inf , +inf] 
     */
    public IloNumVar[] numVarArrayFree(int d1, String name) throws IloException{
        return numVarArray(d1, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, name);
    }
    /**
     * Define a array with free variables: bounds [0 , +inf] 
     */
    public IloNumVar[] numVarArrayPos(int d1, String name) throws IloException{
        return numVarArray(d1, 0, Double.POSITIVE_INFINITY, name);
    }
    /**
     * Define a array with free variables: bounds [-inf , 0] 
     */
    public IloNumVar[] numVarArrayNeg(int d1, String name) throws IloException{
        return numVarArray(d1, Double.NEGATIVE_INFINITY, 0, name);
    }
    /**
     * Define a array with free variables: bounds [-inf , +inf] 
     */
    public IloNumVar[][] numVarArrayFree(int d1, int d2, String name) throws IloException{
        return numVarArray(d1, d2, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, name);
    }
    /**
     * Define a array with free variables: bounds [0 , +inf] 
     */
    public IloNumVar[][] numVarArrayPos(int d1, int d2, String name) throws IloException{
        return numVarArray(d1, d2, 0, Double.POSITIVE_INFINITY, name);
    }
    /**
     * Define a array with free variables: bounds [-inf , 0] 
     */
    public IloNumVar[][] numVarArrayNeg(int d1, int d2,  String name) throws IloException{
        return numVarArray(d1, d2, Double.NEGATIVE_INFINITY, 0, name);
    }
    /**
     * Define a array with free variables: bounds [-inf , +inf] 
     */
    public IloNumVar[][][] numVarArrayFree(int d1, int d2, int d3, String name) throws IloException{
        return numVarArray(d1, d2, d3, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, name);
    }
    /**
     * Define a array with free variables: bounds [0 , +inf] 
     */
    public IloNumVar[][][] numVarArrayPos(int d1, int d2, int d3, String name) throws IloException{
        return numVarArray(d1, d2, d3, 0, Double.POSITIVE_INFINITY, name);
    }
    /**
     * Define a array with free variables: bounds [-inf , 0] 
     */
    public IloNumVar[][][] numVarArrayNeg(int d1, int d2, int d3,  String name) throws IloException{
        return numVarArray(d1, d2, d3, Double.NEGATIVE_INFINITY, 0, name);
    }
    /**
     * Define a array with free variables: bounds [-inf , +inf] 
     */
    public IloNumVar[][][][] numVarArrayFree(int d1, int d2, int d3, int d4, String name) throws IloException{
        return numVarArray(d1, d2, d3, d4, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, name);
    }
    /**
     * Define a array with free variables: bounds [0 , +inf] 
     */
    public IloNumVar[][][][] numVarArrayPos(int d1, int d2, int d3, int d4, String name) throws IloException{
        return numVarArray(d1, d2, d3, d4, 0, Double.POSITIVE_INFINITY, name);
    }
    /**
     * Define a array with free variables: bounds [-inf , 0] 
     */
    public IloNumVar[][][][] numVarArrayNeg(int d1, int d2, int d3, int d4,  String name) throws IloException{
        return numVarArray(d1, d2, d3, d4, Double.NEGATIVE_INFINITY, 0, name);
    }
    /**
     * Define a array with free variables: bounds [-inf , +inf] 
     */
    public IloNumVar[][][][][] numVarArrayFree(int d1, int d2, int d3, int d4, int d5, String name) throws IloException{
        return numVarArray(d1, d2, d3, d4, d5, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, name);
    }
    /**
     * Define a array with free variables: bounds [0 , +inf] 
     */
    public IloNumVar[][][][][] numVarArrayPos(int d1, int d2, int d3, int d4, int d5, String name) throws IloException{
        return numVarArray(d1, d2, d3, d4, d5, 0, Double.POSITIVE_INFINITY, name);
    }
    /**
     * Define a array with free variables: bounds [-inf , 0] 
     */
    public IloNumVar[][][][][] numVarArrayNeg(int d1, int d2, int d3, int d4, int d5, String name) throws IloException{
        return numVarArray(d1, d2, d3, d4, d5, Double.NEGATIVE_INFINITY, 0, name);
    }
    public IloNumVar[] numVarArray(int d1, double min, double max, String name) throws IloException{
        IloNumVar var[] = new IloNumVar[d1];
        for(int i=0; i<d1; i++){
            var[i] = numVar(min, max, name+Index(i, d1));
        }
        return var;
    }
    public IloNumVar[][] numVarArray(int d1, int d2, double min, double max, String name) throws IloException{
        IloNumVar var[][] = new IloNumVar[d1][];
        for(int i=0; i<d1; i++){
            var[i] = numVarArray(d2, min, max, name+Index(i, d1));
        }
        return var;
    }
    public IloNumVar[][][] numVarArray(int d1, int d2, int d3, double min, double max, String name) throws IloException{
        IloNumVar var[][][] = new IloNumVar[d1][][];
        for(int i=0; i<d1; i++){
            var[i] = numVarArray(d2, d3, min, max, name+Index(i, d1));
        }
        return var;
    }
    public IloNumVar[][][][] numVarArray(int d1, int d2, int d3, int d4, double min, double max, String name) throws IloException{
        IloNumVar var[][][][] = new IloNumVar[d1][][][];
        for(int i=0; i<d1; i++){
            var[i] = numVarArray(d2, d3, d4, min, max, name+Index(i, d1));
        }
        return var;
    }
    public IloNumVar[][][][][] numVarArray(int d1, int d2, int d3, int d4, int d5, double min, double max, String name) throws IloException{
        IloNumVar var[][][][][] = new IloNumVar[d1][][][][];
        for(int i=0; i<d1; i++){
            var[i] = numVarArray(d2, d3, d4, d5, min, max, name+Index(i, d1));
        }
        return var;
    }
    public IloIntVar[] intVarArray(int d1, int min, int max, String name) throws IloException{
        IloIntVar var[] = new IloIntVar[d1];
        for(int i=0; i<d1; i++){
            var[i] = intVar(min, max, name+Index(i, d1));
        }
        return var;
    }
    public IloIntVar[][] intVarArray(int d1, int d2, int min, int max, String name) throws IloException{
        IloIntVar var[][] = new IloIntVar[d1][];
        for(int i=0; i<d1; i++){
            var[i] = intVarArray(d2, min, max, name+Index(i, d1));
        }
        return var;
    }
    public IloIntVar[][][] intVarArray(int d1, int d2, int d3, int min, int max, String name) throws IloException{
        IloIntVar var[][][] = new IloIntVar[d1][][];
        for(int i=0; i<d1; i++){
            var[i] = intVarArray(d2, d3, min, max, name+Index(i, d1));
        }
        return var;
    }
    public IloIntVar[][][][] intVarArray(int d1, int d2, int d3, int d4, int min, int max, String name) throws IloException{
        IloIntVar var[][][][] = new IloIntVar[d1][][][];
        for(int i=0; i<d1; i++){
            var[i] = intVarArray(d2, d3, d4, min, max, name+Index(i, d1));
        }
        return var;
    }
    public IloIntVar[] boolVarArray(int d1, String name) throws IloException{
        IloIntVar var[] = new IloIntVar[d1];
        for(int i=0; i<d1; i++){
            var[i] = boolVar(name+Index(i, d1));
        }
        return var;
    }
    public IloIntVar[][] boolVarArray(int d1, int d2, String name) throws IloException{
        IloIntVar var[][] = new IloIntVar[d1][];
        for(int i=0; i<d1; i++){
            var[i] = boolVarArray(d2, name+Index(i, d1));
        }
        return var;
    }
    public IloIntVar[][][] boolVarArray(int d1, int d2, int d3, String name) throws IloException{
        IloIntVar var[][][] = new IloIntVar[d1][][];
        for(int i=0; i<d1; i++){
            var[i] = boolVarArray(d2, d3, name+Index(i, d1));
        }
        return var;
    }
    public IloIntVar[][][][] boolVarArray(int d1, int d2, int d3, int d4, String name) throws IloException{
        IloIntVar var[][][][] = new IloIntVar[d1][][][];
        for(int i=0; i<d1; i++){
            var[i] = boolVarArray(d2, d3, d4,name+Index(i, d1));
        }
        return var;
    }
    public IloIntVar[][][][][] boolVarArray(int d1, int d2, int d3, int d4, int d5, String name) throws IloException{
        IloIntVar var[][][][][] = new IloIntVar[d1][][][][];
        for(int i=0; i<d1; i++){
            var[i] = boolVarArray(d2, d3, d4, d5, name+Index(i, d1));
        }
        return var; 
    }

    public IloNumVar[][] numVarArray(int d1, int d2, double min, double max) throws IloException{
        IloNumVar var[][] = new IloNumVar[d1][];
        for(int i=0; i<d1; i++){
            var[i] = numVarArray(d2, min, max);
        }
        return var;
    }
    public IloNumVar[][][] numVarArray(int d1, int d2, int d3, double min, double max) throws IloException{
        IloNumVar var[][][] = new IloNumVar[d1][][];
        for(int i=0; i<d1; i++){
            var[i] = numVarArray(d2, d3, min, max);
        }
        return var;
    }
    public IloNumVar[][][][] numVarArray(int d1, int d2, int d3, int d4,double min, double max) throws IloException{
        IloNumVar var[][][][] = new IloNumVar[d1][][][];
        for(int i=0; i<d1; i++){
            var[i] = numVarArray(d2, d3, d4, min, max);
        }
        return var;
    }
    public IloNumVar[][][][][] numVarArray(int d1, int d2, int d3, int d4, int d5, double min, double max) throws IloException{
        IloNumVar var[][][][][] = new IloNumVar[d1][][][][];
        for(int i=0; i<d1; i++){
            var[i] = numVarArray(d2, d3, d4, d5, min, max);
        }
        return var;
    }
    
    @Override
    public double getValue(IloNumExpr V, boolean ignoreUnknown) throws IloException{
        if(V==null){return Double.NaN;}
        
        if(ignoreUnknown){
            try{
                return super.getValue(V);
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
    public double[] getValues(IloNumExpr V[]) throws UnknownObjectException, IloException{
        if(V==null){return null;}
        
        return getValues(V, false);
    }
    public double[] getValues(IloNumExpr V[], boolean ignoreUnknown) throws UnknownObjectException, IloException{
        if(V==null){return null;}
        
        double X[] = new double[V.length];
        for(int i=0; i<V.length; i++){
            X[i] = getValue(V[i], ignoreUnknown);
        }
        return X;
    }
    public double[][] getValues(IloNumExpr V[][]) throws UnknownObjectException, IloException{
        if(V==null){return null;}
        
        return getValues(V, false);
    }
    public double[][] getValues(IloNumExpr V[][], boolean ignoreUnknown) throws UnknownObjectException, IloException{
        if(V==null){return null;}
        
        double X[][] = new double[V.length][];
        for(int i=0; i<V.length; i++){
            X[i] = getValues(V[i], ignoreUnknown);
        }
        return X;
    }
    public double[][][] getValues(IloNumExpr V[][][]) throws UnknownObjectException, IloException{
        if(V==null){return null;}
        
        return getValues(V, false);
    }
    public double[][][] getValues(IloNumExpr V[][][], boolean ignoreUnknown) throws UnknownObjectException, IloException{
        if(V==null){return null;}
        
        double X[][][] = new double[V.length][][];
        for(int i=0; i<V.length; i++){
            X[i] = getValues(V[i], ignoreUnknown);
        }
        return X;
    }
    public double[][][][] getValues(IloNumExpr V[][][][]) throws UnknownObjectException, IloException{
        if(V==null){return null;}
        
        return getValues(V, false);
    }
    public double[][][][] getValues(IloNumExpr V[][][][], boolean ignoreUnknown) throws UnknownObjectException, IloException{
        if(V==null){return null;}
        
        double X[][][][] = new double[V.length][][][];
        for(int i=0; i<V.length; i++){
            X[i] = getValues(V[i], ignoreUnknown);
        }
        return X;
    }

    public IloNumExpr Sum(IloNumExpr sum, IloNumExpr exp, double b) throws IloException{
        if(sum==null){
            return sum(exp,b);
        }else if(exp==null){
            return sum(sum,b);
        }else {
            return sum(sum(sum, exp),b);
        }
    }
    public IloNumExpr Sum(IloNumExpr sum, IloNumExpr exp) throws IloException{
        if(sum==null){
            return exp;
        }else if(exp==null){
            return sum;
        }else{
            return sum(sum, exp);
        }
    }
    public IloNumExpr Sum(IloNumExpr... exp) throws IloException{
        IloNumExpr sum = null;
        for(IloNumExpr e : exp){
            sum = Sum(sum, e);
        }
        return sum;
    }
//    public IloNumExpr Sum(IloNumExpr M[]) throws IloException{
//        return sum(M);
//    }
    public IloNumExpr Sum(IloNumExpr M[][]) throws IloException{
        IloNumExpr aux[] = new IloNumExpr[M.length];
        for(int i=0; i<M.length; i++){
            aux[i] = Sum(M[i]);
        }
        return sum(aux);
    }
    public IloNumExpr Sum(IloNumExpr M[][][]) throws IloException{
        IloNumExpr aux[] = new IloNumExpr[M.length];
        for(int i=0; i<M.length; i++){
            aux[i] = Sum(M[i]);
        }
        return sum(aux);
    }
    public IloNumExpr Sum(IloNumExpr M[][][][]) throws IloException{
        IloNumExpr aux[] = new IloNumExpr[M.length];
        for(int i=0; i<M.length; i++){
            aux[i] = Sum(M[i]);
        }
        return sum(aux);
    }

    public void addSubject(Object ...obj) throws IloException{
        int signal = 1;
        int opRel = -1;
        LinkedList<IloNumExpr> A = new LinkedList<IloNumExpr>();
        double sum = 0;

        for(Object o:obj){
            if(o instanceof String){
                String s = (String)o;
                if(s.equals("+")){
                    signal = 1;
                }else if(s.equals("-")){
                    signal = -1;
                }else if(s.equals("Le") && opRel==-1){
                    signal = 1;
                    opRel = 0;
                }else if(s.equals("Eq") && opRel==-1){
                    signal = 1;
                    opRel = 1;
                }else if(s.equals("Ge") && opRel==-1){
                    signal = 1;
                    opRel = 2;
                }else{
                    throw new IloException("String s = "+s+" not valid, opRel="+opRel);
                }
            }else if(o instanceof IloNumExpr[]){
                IloNumExpr[] var = (IloNumExpr[])o;
                if(opRel==-1){
                    A.addLast(prod(signal, sum(var)));
                }else{
                    A.addLast(prod( -signal , sum(var)));
                }
                signal = 1;
            }else if(o instanceof IloNumExpr){
                IloNumExpr var = (IloNumExpr)o;
                if(opRel==-1){
                    A.addLast(prod(signal, var));
                }else{
                    A.addLast(prod( -signal , var));
                }
                signal = 1;
            }else if(o instanceof Double){
                Double constante = (Double)o;
                if(opRel==-1){
                    sum -= signal*constante;
                }else{
                    sum += signal*constante;
                }
                signal = 1;
            }else if(o instanceof Integer){
                Integer constante = (Integer)o;
                if(opRel==-1){
                    sum -= signal*constante;
                }else{
                    sum += signal*constante;
                }
                signal = 1;
            }else{
                throw new IloException("Object o = "+o+" not known");
            }
        }
        if(opRel==0){
            addLe(sum(A.toArray(new IloNumExpr[A.size()])), sum);
        }else if(opRel==1){
            addEq(sum(A.toArray(new IloNumExpr[A.size()])), sum);
        }else if(opRel==2){
            addGe(sum(A.toArray(new IloNumExpr[A.size()])), sum);
        }else{
            throw new IloException("Operator relatinal not found");
        }
    }
    
    
    public final IloNumExpr Not(IloNumExpr P) throws IloException {
        return sum(1, prod(-1, P));
    }
//    public final IloNumVar Not(String name, IloNumVar P) throws IloException {
//        IloNumVar y = numVar(0.0, 1.0, name);
//        if(name==null){
//            addEq(y, sum(1, prod(-1, P)));
//        }else{
//            addEq(y, sum(1, prod(-1, P)), name+":not");
//        }
//        return y;
//    }
    
    public final IloNumVar And(IloNumExpr P1, IloNumExpr P2) throws IloException {
        return And(null, P1, P2);
    }
    public final IloNumVar Or(IloNumExpr P1, IloNumExpr P2) throws IloException {
        return Or(null, P1, P2);
    }
    public final IloNumVar XOr(IloNumExpr P1, IloNumExpr P2) throws IloException {
        return XOr(null, P1, P2);
    }
    public final IloNumVar IF_Then(IloNumExpr P1, IloNumExpr P2) throws IloException {
        return IF_Then(null, P1, P2);
    }
    public final IloNumVar IF_Only(IloNumExpr P1, IloNumExpr P2) throws IloException {
        return IF_Only(null, P1, P2);
    }
    public final IloNumVar And(String name, IloNumExpr P1, IloNumExpr P2) throws IloException {
        IloNumVar y = numVar(0.0, 1.0, name);
        addAnd(name, y, P1, P2);
        return y;
    }
    public final IloNumVar Or(String name, IloNumExpr P1, IloNumExpr P2) throws IloException {
        IloNumVar y = numVar(0.0, 1.0, name);
        addOr(name, y, P1, P2);
        return y;
    }
    public final IloNumVar XOr(String name, IloNumExpr P1, IloNumExpr P2) throws IloException {
        IloNumVar y = numVar(0.0, 1.0, name);
        addXOr(name, y, P1, P2);
        return y;
    }
    public final IloNumVar IF_Then(String name, IloNumExpr P1, IloNumExpr P2) throws IloException {
        IloNumVar y = numVar(0.0, 1.0, name);
        addIF_Then(name, y, P1, P2);
        return y;
    }
    public final IloNumVar IF_Only(String name, IloNumExpr P1, IloNumExpr P2) throws IloException {
        IloNumVar y = numVar(0.0, 1.0, name);
        addIF_Only(name, y, P1, P2);
        return y;
    }
    public final void addAnd(String name, IloNumExpr y, IloNumExpr P1, IloNumExpr P2) throws IloException {
//        add_00_0(name, y, P1, P2);
//        add_01_0(name, y, P1, P2);
//        add_10_0(name, y, P1, P2);
//        add_11_1(name, y, P1, P2);
        if(name!=null){
            addGe(y, sumArg(+1,P1, +1,P2, -1), name+":00_1");
            addLe(y, P1,                    name+":10_1");
            addLe(y, P2,                    name+":01_1");
        }else{
            addGe(y, sumArg(+1,P1, +1,P2, -1));
            addLe(y, P1);
            addLe(y, P2);
        }
    }
    public final IloNumVar And(IloNumExpr... Pi) throws IloException {
        return And(null, Pi);
    }
    public final IloNumVar And(String name, IloNumExpr... Pi) throws IloException {
        IloNumVar y = numVar(0.0, 1.0, name);
        addAnd(name, y, Pi);
        return y;
    }
    public final void addAnd(String name, IloNumExpr y, IloNumExpr... Pi) throws IloException {
//        IloNumVar wi[] = numVarArray(Pi.length, 0.0, 1.0, "w");
//        addEq(wi[0], Pi[0], name+",w[0].p[0]");
//        for(int i=1; i<Pi.length; i++){
//            And(name+",w["+i+"].w["+(i-1)+"]&p["+i+"]", wi[i], wi[i-1], Pi[i]);
//        }
        //addEq(y, wi[Pi.length-1], "y.w["+(Pi.length-1)+"]");
        addGe(y, sum(1-Pi.length, sum(Pi)), name+":Tand_1");
        for(int i=0; i<Pi.length; i++){
            addLe(y, Pi[i], name+":Tand_0");
        }
    }
    public final void addOr(String name, IloNumExpr y, IloNumExpr P1, IloNumExpr P2) throws IloException {
//        add_00_0(name, y, P1, P2);
//        add_01_1(name, y, P1, P2);
//        add_10_1(name, y, P1, P2);
//        add_11_1(name, y, P1, P2);
        if(name!=null){
            addLe(y, sumArg(+1,P1, +1,P2),  name+":00_1");
            addGe(y, P1,                    name+":10_1");
            addGe(y, P2,                    name+":01_1");
        }else{
            addLe(y, sumArg(+1,P1, +1,P2));
            addGe(y, P1);
            addGe(y, P2);
        }
    }
    public final IloNumVar Or(IloNumExpr... Pi) throws IloException {
        return Or(null, Pi);
    }
    public final IloNumVar Or(String name, IloNumExpr... Pi) throws IloException {
        IloNumVar y = numVar(0.0, 1.0, name);
        addOr(name, y, Pi);
        return y;
    }
    public final void addOr(String name, IloNumVar y, IloNumExpr... Pi) throws IloException {
//        IloNumVar wi[] = numVarArray(Pi.length, 0.0, 1.0, "w");
//        addEq(wi[0], Pi[0], name+",w[0].p[0]");
//        for(int i=1; i<Pi.length; i++){
//            And(name+",w["+i+"].w["+(i-1)+"]&p["+i+"]", wi[i], wi[i-1], Pi[i]);
//        }
        //addEq(y, wi[Pi.length-1], "y.w["+(Pi.length-1)+"]");
        addLe(y, sum(Pi), name+":Tor_1");
        for(int i=0; i<Pi.length; i++){
            addGe(y, Pi[i], name+":Tor_0");
        }
    }
    
    public final void addXOr(String name, IloNumExpr y, IloNumExpr P1, IloNumExpr P2) throws IloException {
        add_00_0(name, y, P1, P2);
        add_01_1(name, y, P1, P2);
        add_10_1(name, y, P1, P2);
        add_11_0(name, y, P1, P2);
    }
    public final void addIF_Then(String name, IloNumExpr y, IloNumExpr P1, IloNumExpr P2) throws IloException {
        add_00_1(name, y, P1, P2);
        add_01_1(name, y, P1, P2);
        add_10_0(name, y, P1, P2);
        add_11_1(name, y, P1, P2);
    }
    public final void addIF_Only(String name, IloNumExpr y, IloNumExpr P1, IloNumExpr P2) throws IloException {
        add_00_1(name, y, P1, P2);
        add_01_0(name, y, P1, P2);
        add_10_0(name, y, P1, P2);
        add_11_1(name, y, P1, P2);
    }
    
    private final void add_00_0(String name, IloNumExpr y, IloNumExpr p1, IloNumExpr p2) throws IloException {
        if(name!=null){
            addLe(y, sum(p1, p2), name+":00_0");
        }else{
            addLe(y, sum(p1, p2));
        }
    }
    private final void add_01_0(String name, IloNumExpr y, IloNumExpr p1, IloNumExpr p2) throws IloException {
        if(name!=null){
            addLe(y, sumArg(p1, -1,p2, +1), name+":01_0");
        }else{
            addLe(y, sumArg(p1, -1,p2, +1));
        }
    }
    private final void add_10_0(String name, IloNumExpr y, IloNumExpr p1, IloNumExpr p2) throws IloException {
        if(name!=null){
            addLe(y, sumArg(p2, -1,p1, +1), name+":10_0");
        }else{
            addLe(y, sumArg(p2, -1,p1, +1));
        }
    }
    private final void add_11_0(String name, IloNumExpr y, IloNumExpr p1, IloNumExpr p2) throws IloException {
        if(name!=null){
            addLe(y, sumArg(-1,p1, -1,p2, +2), name+":11_0");
        }else{
            addLe(y, sumArg(-1,p1, -1,p2, +2));
        }
    }
    private final void add_00_1(String name, IloNumExpr y, IloNumExpr p1, IloNumExpr p2) throws IloException {
        if(name!=null){
            addGe(y, sumArg(-1,p1, -1,p2, +1), name+":00_1");
        }else{
            addGe(y, sumArg(-1,p1, -1,p2, +1));
        }
    }
    private final void add_01_1(String name, IloNumExpr y, IloNumExpr p1, IloNumExpr p2) throws IloException {
        if(name!=null){
            addGe(y, sumArg(-1,p1, +1,p2), name+":01_1");
        }else{
            addGe(y, sumArg(-1,p1, +1,p2));
        }
    }
    private final void add_10_1(String name, IloNumExpr y, IloNumExpr p1, IloNumExpr p2) throws IloException {
        if(name!=null){
            addGe(y, sumArg(+1,p1, -1,p2), name+":10_1");
        }else{
            addGe(y, sumArg(+1,p1, -1,p2));
        }
    }
    private final void add_11_1(String name, IloNumExpr y, IloNumExpr p1, IloNumExpr p2) throws IloException {
        if(name!=null){
            addGe(y, sumArg(+1,p1, +1,p2, -1), name+":00_1");
        }else{
            addGe(y, sumArg(+1,p1, +1,p2, -1));
        }
    } 
    
    public final IloNumExpr sumArg(Object ...args) throws IloException{
//        System.out.println("-------------------------------------------------"); 
//        for(int i=0; i<args.length; i++){
//            System.out.printf("%s,", args[i]);
//        }
//        System.out.println(); 
        IloNumExpr exp = null;
        for(int i=0; i<args.length; i++){
            if(args[i] instanceof Double){
                double coef = (Double)args[i];
                if(i+1<args.length){
                    IloNumExpr var = (IloNumExpr)args[i+1];
                    i++;
                    if(exp==null){
                        exp = prod(coef, var);
                    }else{
                        exp = sum(prod(coef, var), exp);
                    }
                }else{
                    exp = sum(exp, coef);
                }
            }else if(args[i] instanceof Integer){
                double coef = (Integer)args[i];
                if(i+1<args.length){
                    IloNumExpr var = (IloNumExpr)args[i+1];
                    i++;
                    if(exp==null){
                        exp = prod(coef, var);
                    }else{
                        exp = sum(prod(coef, var), exp);
                    }
                }else{
                    exp = sum(exp, coef);
                }
            }else if(args[i] instanceof IloNumExpr){
                IloNumExpr var = (IloNumExpr)args[i];
                if(exp==null){
                    exp = var;
                }else{
                    exp = sum(var, exp);
                }
            }else{
                throw new IloException("arg["+i+"] = '"+args[i]+"' is not valid");
            }
            //System.out.printf("exp ~ %s\n", exp);
            
        }
        
        //System.out.println("exp = "+exp);
        return exp;
    }
    
    /**
     * <b>If</b> y = 1 <b>then</b> Zi = 1 for all (i);
     * @param name
     * @param y
     * @param Zi
     * @throws IloException 
     */
    public final void addIF_Y_Them_Zi(String name, IloNumExpr y, IloNumExpr ...Zi) throws IloException {
        for(IloNumExpr z : Zi){
            addLe(y, z, name);    
        }
    }
    
    /**
     * <b>If</b> y = 1 <b>then</b> constraints &le 0;
     * @param name
     * @param y
     * @param constraints
     * @throws IloException 
     */
    public final void addIF_Y_Them_Le(String name, IloNumExpr y, IloNumExpr ...constraints) throws IloException {
        for(IloNumExpr exp : constraints){
            //addLe(exp, sum(bigM-epsilon, prod(-bigM, y)), name);    //Ax - b <= M(1-y) - e  || (M-e) - My
            addLe(exp, sum(bigM, prod(-bigM, y)), name);    //Ax - b <= M(1-y)  || M - My
        }
    }
    /**
     * <b>If</b> y = 1 <b>then</b> constraints &ge 0;
     * @param name
     * @param y
     * @param constraints
     * @throws IloException 
     */
    public final void addIF_Y_Them_Ge(String name, IloNumExpr y, IloNumExpr ...constraints) throws IloException {
        for(IloNumExpr exp : constraints){
            //addGe(exp, sum(epsilon-bigM, prod(+bigM, y)), name);    //Ax - b >= M(y-1) + e  || (e-M) + My
            addGe(exp, sum(-bigM, prod(+bigM, y)), name);    //Ax - b >= M(y-1)  || -M + My
        }
    }
    /**
     * <b>If</b> y = 1 <b>then</b> constraints = 0;
     * @param name
     * @param y
     * @param constraints
     * @throws IloException 
     */
    public final void addIF_Y_Them_Eq(String name, IloNumExpr y, IloNumExpr ...constraints) throws IloException {
        //Ax - b <= M(1-y)+e
        //b - Ax <= M(1-y)+e
        for(IloNumExpr exp : constraints){
            addLe(exp, sum(+bigM+epsilon, prod(-bigM, y)), name);    //Ax - b <= M(1-y) + e  || +(M+e) - My
            addGe(exp, sum(-bigM-epsilon, prod(+bigM, y)), name);    //Ax - b >= M(y-1) - e  || -(M+e) + My
        }
    }
    
    /**
     * <b>If</b> y = 1 <b>then</b> <i>constraints</i> = 0;
     * @param name
     * @param LB is the minimum possible value for <i>constraints</i>
     * @param UB is the maximum possible value for <i>constraints</i>
     * @param epsilon a small possitive value, sujestion 1e-4
     * @param y is a binary proposition
     * @param constraints is set of expression
     * @throws IloException 
     */
    public final void addIF_Y_Them_Eq(String name, final double LB, final double UB, final double epsilon, IloNumExpr y, IloNumExpr ...constraints) throws IloException {
        //Ax - b <= UB*(1-y)+e
        //Ax - b >= LB*(1-y)-e
        for(IloNumExpr exp : constraints){
            addLe(exp, sum(+UB+epsilon, prod(-UB, y)), name);    //Ax - b <= UB*(1-y) + e  || +(UB+e) - UB*y
            addGe(exp, sum(+LB-epsilon, prod(-LB, y)), name);    //Ax - b >= LB*(1-y) - e  || +(LB-e) - LB*y
        }
    }
    
    /**
     * <b>If</b> y = 1 <b>then</b> exp1 &le q &le exp2;
     * @param name
     * @param y
     * @param exp1
     * @param q
     * @param exp2
     * @throws IloException 
     */
    public void addIF_Y_Them_Between(String name, IloNumExpr y, IloNumExpr exp1, double q, IloNumExpr exp2) throws IloException {
        addIF_Y_Them_Le(name, y, sum(exp1, -q));
        addIF_Y_Them_Ge(name, y, sum(exp2, -q));
    }
    /**
     * <b>If</b> y = 1 <b>then</b> lb &le exp &le ub;
     * <br><br><b>Complexity (B,C,R):</b>  (-, -, 2)    <br>
     * <pre>
     * This call add two new restrictions to the model:
     * 
     *      exp &le ub + M(1-y)
     *      exp &ge lb + M(y-1)
     * 
     * where M is a big positive value
     * </pre>
     *
     * @param name
     * @param y
     * @param lb
     * @param exp
     * @param ub
     * @throws IloException 
     */
    public void addIF_Y_Them_Between(String name, IloNumExpr y, double lb, IloNumExpr exp, double ub) throws IloException {
        addIF_Y_Them_Le(name, y, sum(exp, -ub));
        addIF_Y_Them_Ge(name, y, sum(exp, -lb));
    }
    
    /**
     * <b>If</b> y = 1 <b>then</b> z = w;
     * @param name
     * @param y is a binary proposition
     * @param z is a binary proposition
     * @param w is a binary proposition
     * @throws IloException 
     */
    public final void addIF_Y_Them_Z_Eq_W(String name, IloNumExpr y, IloNumExpr z, IloNumExpr w) throws IloException {
        addIF_Y_Them_Eq(name, -1, +1, epsilon, y, sum(z, prod(-1, w)));
    }
    /**
     * <b>If</b> y = 1 <b>then</b> z = value;
     * @param name
     * @param y is a proposition variable
     * @param z is a proposition variable
     * @param value is a parameter value {true or false}
     * @throws IloException 
     */
    public final void addIF_Y_Them_Z_Eq_W(String name, IloNumExpr y, IloNumExpr z, final boolean value) throws IloException {
        if(value){
            addGe(z, y, name);
        }else{
            addLe(sum(z, y), 1, name);
        }
    }
    
    public final void addIF_Y_Them_StayIn(String name, IloNumExpr y, Hyperplane hyperplans[], double risks[], IloNumVar x[]) throws IloException {
        for(int i=0; i<hyperplans.length; i++){
            //final double risk = approach.unc.risk_allocation(t, fixed_delta);
            //final double risk = system.unc.risk_allocation(t, fixed_delta[c], obs.hyperplans[i].a);
            //c_(i,t) (δ_it) ≤ a_i^T*x_t - b_i 

            IloNumExpr exp = hyperplans[i].scalProd(this, x);
            exp = sum(exp, -hyperplans[i].b - (risks==null?0:risks[i]));
            addIF_Y_Them_Ge(name, y, exp);
        }
    }
    public final void addIF_Y_Them_StayOut(String name, IloNumExpr y, Hyperplane hyperplans[], double risks[], IloNumVar x[]) throws IloException {
        IloNumVar W[] = boolVarArray(hyperplans.length, name+".w");
        for(int i=0; i<hyperplans.length; i++){
            //final double risk = approach.unc.risk_allocation(t, fixed_delta);
            //final double risk = system.unc.risk_allocation(t, fixed_delta[c], obs.hyperplans[i].a);
            //c_(i,t) (δ_it) ≤ a_i^T*x_t - b_i 
            
            IloNumExpr exp = hyperplans[i].scalProd(this, x);
            exp = sum(exp, -hyperplans[i].b - (risks==null?0:risks[i]));
            addIF_Y_Them_Ge(name, W[i], exp);
        }
        addGe(sum(W), y, "(4.91).sum_ge_y");
    }
 
    /**
     * sum{ Z(t) and Z(t-1) } = 1
     * @param name
     * @param y
     * @param hyperplans
     * @param risks
     * @param x
     * @throws IloException 
     */
    public final void addIF_Y_Them_StayOut(String name, IloNumExpr y[], Hyperplane hyperplans[], double risks[], IloNumVar x[][]) throws IloException {
        IloNumVar W[][] = boolVarArray(y.length, hyperplans.length, name+".w");
        for(int t=1; t<y.length; t++){
            IloNumExpr sum = null;
            for (int i = 0; i < hyperplans.length; i++) {
                sum = SumProd(sum, 1.0, And(name+".And", W[t][i], W[t-1][i]));
            }
            addGe(sum, y[t], name+".Or");
        }
        for(int t=0; t<y.length; t++){
            for(int i=0; i<hyperplans.length; i++){
                IloNumExpr exp = hyperplans[i].scalProd(this, x[t]);
                exp = sum(exp, -hyperplans[i].b - (risks==null?0:risks[i]));
                addIF_Y_Them_Ge(name, W[t][i], exp);
            }
        }
    }
    public final void addIF_Y_Them_StayOut(String name, IloNumExpr y[], Hyperplane hyperplans[], double risks[], IloNumVar x[][], IloNumExpr erf_inv[], final int W_length) throws IloException {
        System.out.println("y = "+y.length);
        System.out.println("w = "+W_length);
        System.out.println("x = "+x.length);
        
        
        //x>=w>=y
        IloNumVar W[][] = boolVarArray(W_length+1, hyperplans.length, name+".w");
        for(int t=1; t<W.length; t++){
            IloNumExpr sum = null;
            for (int i = 0; i < hyperplans.length; i++) {
                sum = SumProd(sum, 1.0, And(name+".And", W[t][i], W[t-1][i]));
            }
            
            int s = (t*(y.length-1))/(W.length-1);
            
            addGe(sum, y[s], name+".Or");
            System.out.printf("t = %d -> %d\n", t, s);
        }
//        for(int t=0; t<W.length; t++){
//            addGe(sum(W[t]), y[(t*(y.length-1))/(W.length-1)], "(4.91).sum_ge_y");
//            System.out.printf("t = %d -> %d\n", t, ((t*(y.length-1))/(W.length-1)));
//        }
        
        for(int t=0; t<x.length; t++){
            int s = (t*(W.length-1))/(x.length-1);
            for(int i=0; i<hyperplans.length; i++){
                IloNumExpr exp = hyperplans[i].scalProd(this, x[t]);
                if(risks==null && Math.abs(risks[i])>1e-6){
                    exp = sum(exp, prod(-risks[i], erf_inv[t]));
                }
                exp = sum(exp, -hyperplans[i].b);
                addIF_Y_Them_Ge(name, W[s][i], exp);
                //addIF_Y_Them_Ge(name, W[t][i], exp);
            }
            System.out.printf("t = %d -> %d\n", t, s);
        }
    }
    
    public final void addIF_Y_Them_StayOut(String name, IloNumExpr y[], Hyperplane hyperplans[], double sigma[], IloNumVar x[][], final int W_length) throws IloException {
        System.out.println("y = "+y.length);
        System.out.println("w = "+W_length);
        System.out.println("x = "+x.length);
        
        
        //x>=w>=y
        IloNumVar W[][] = boolVarArray(W_length+1, hyperplans.length, name+".w");
        //int twy[] = new int[W_length];
        for(int t=1; t<W.length; t++){
            IloNumExpr sum = null;
            for (int i = 0; i < hyperplans.length; i++) {
                sum = SumProd(sum, 1.0, And(name+".And", W[t][i], W[t-1][i]));
            }
            
            int s = (t*(y.length-1))/(W.length-1);
            
            addGe(sum, y[s], name+".Or");
            System.out.printf("t = %d -> %d\n", t, s);
        }
//        for(int t=0; t<W.length; t++){
//            addGe(sum(W[t]), y[(t*(y.length-1))/(W.length-1)], "(4.91).sum_ge_y");
//            System.out.printf("t = %d -> %d\n", t, ((t*(y.length-1))/(W.length-1)));
//        }
        
        for(int t=0; t<x.length; t++){
            int s = (t*(W.length-1))/(x.length-1);
            for(int i=0; i<hyperplans.length; i++){
                IloNumExpr exp = hyperplans[i].scalProd(this, x[t]);
                exp = sum(exp, -hyperplans[i].b - (sigma==null?0:sigma[i]));
                addIF_Y_Them_Ge(name, W[s][i], exp);
                //addIF_Y_Them_Ge(name, W[t][i], exp);
            }
            System.out.printf("t = %d -> %d\n", t, s);
        }
    }
    
    /**
     * <b>Attention, requires the creation of new binary variables for each constraint</b><br>
     * <b>If</b> constraints &le 0 <b>then</b> y = 1;
     * @param name
     * @param y
     * @param constraints
     * @throws IloException 
     */
    public final void addIF_Le_Them_Y(String name, IloNumExpr y, IloNumExpr ...constraints) throws IloException {
        IloIntVar zi[] = boolVarArray(constraints.length, "z");
        for(int i=0; i<constraints.length; i++){
            addIF_Le_Them_Y(name+"["+i+"]", zi[i], constraints[i]);
        }
        addAnd(name+".and", y, zi);
    }
    /**
     * <b>If</b> constraint &le 0 <b>then</b> y = 1;
     * @param name
     * @param y
     * @param constraint
     * @throws IloException 
     */
    public final void addIF_Le_Them_Y(String name, IloNumExpr y, IloNumExpr constraint) throws IloException {
        addGe(constraint, sum(+epsilon, prod(-bigM, y)), name);    //Ax - b >= e - My
    }
    /**
     * <b>If</b> constraint &lt 0 <b>then</b> y = 1;
     * @param name
     * @param y
     * @param constraint
     * @throws IloException 
     */
    public final void addIF_Lt_Them_Y(String name, IloNumExpr y, IloNumExpr constraint) throws IloException {       
        addGe(constraint, prod(-bigM, y), name);    //Ax < b --> y     |    Ax - b >= - My
    }
    /**
     * <b>If</b> constraint &ge 0 <b>then</b> y = 1;
     * @param name
     * @param y
     * @param constraint
     * @throws IloException 
     */
    public final void addIF_Ge_Them_Y(String name, IloNumExpr y, IloNumExpr constraint) throws IloException {
        addLe(constraint, sum(-epsilon, prod(+bigM, y)), name);    //Ax - b <= -e + My
    }
    /**
     * <b>If</b> constraint &gt 0 <b>then</b> y = 1;
     * @param name
     * @param y
     * @param constraint
     * @throws IloException 
     */
    public final void addIF_Gt_Them_Y(String name, IloNumExpr y, IloNumExpr constraint) throws IloException {
        addLe(constraint, prod(+bigM, y), name);    //Ax > b --> y     |    Ax - b <= My
    }

    /**
     * <b>Attention, requires the creation of new binary variables for each constraint</b><br>
     * <b>If</b> constraints &le 0 <b>then</b> y = 1;
     * @param name
     * @param y
     * @param constraints
     * @throws IloException 
     */
    public final void addIF_Eq_Them_Y(String name, IloNumExpr y, IloNumExpr ...constraints) throws IloException {
        IloNumExpr exp[] = new IloNumExpr[2*constraints.length];
        for(int i=0; i<constraints.length; i++){
            exp[2*i] = constraints[i];
            exp[2*i+1] = prod(-1,constraints[i]);
        }
        addIF_Le_Them_Y(name, y, exp);
    }
    public final void addIF_Eq_Them_Y(String name, IloNumExpr y, IloNumExpr exp) throws IloException {
        IloIntVar z = boolVar("z");
        IloIntVar w = boolVar("w");
        addIF_Le_Them_Y(name, z, exp);
        addIF_Ge_Them_Y(name, w, exp);
        addAnd(name+".and", y, z, w);
    }

    /**
     * <b>If</b> exp1 &le q &le exp2 <b>then</b> y = 1;
     * @param name
     * @param exp1
     * @param q
     * @param exp2
     * @return y
     * @throws IloException 
     */
    public IloNumVar IF_Between_Them_Y(String name, IloNumExpr exp1, double q, IloNumExpr exp2) throws IloException {
        IloNumVar y = numVar(0, 1, name+".y");
        addIF_Between_Them_Y(name, y, exp1, q, exp2);
        return y;
    }
    /**
     * <b>If</b> exp1 &le q &le exp2 <b>then</b> y = 1;
     * @param name
     * @param y
     * @param exp1
     * @param q
     * @param exp2
     * @throws IloException 
     */
    public void addIF_Between_Them_Y(String name, IloNumExpr y, IloNumExpr exp1, double q, IloNumExpr exp2) throws IloException {
        IloIntVar z = boolVar("z");
        IloIntVar w = boolVar("w");
        addIF_Le_Them_Y(name, z, sum(exp1, -q));
        addIF_Ge_Them_Y(name, w, sum(exp2, -q));
        addAnd(name+".and", y, z, w);
    }
    /**
     * <b>If</b> exp1 &le q &lt exp2 <b>then</b> y = 1;
     * @param name
     * @param exp1
     * @param q
     * @param exp2
     * @return y
     * @throws IloException 
     */
    public IloNumVar IF_BetweenLT_Them_Y(String name, IloNumExpr exp1, double q, IloNumExpr exp2) throws IloException {
        IloNumVar y = numVar(0, 1, name+".y");
        addIF_BetweenLT_Them_Y(name, y, exp1, q, exp2);
        return y;
    }
    /**
     * <b>If</b> exp1 &le q &lt exp2 <b>then</b> y = 1;
     * @param name
     * @param y
     * @param exp1
     * @param q
     * @param exp2
     * @throws IloException 
     */
    public void addIF_BetweenLT_Them_Y(String name, IloNumExpr y, IloNumExpr exp1, double q, IloNumExpr exp2) throws IloException {
        IloIntVar z = boolVar("z");
        IloIntVar w = boolVar("w");
        addIF_Le_Them_Y(name, z, sum(exp1, -q));
        addIF_Ge_Them_Y(name, w, sum(exp2, -q-epsilon));
        addAnd(name+".and", y, z, w);
    }
    /**
     * <b>If</b> lb &le exp &le ub <b>then</b> y = 1;
     * @param name
     * @param lb
     * @param exp
     * @param ub
     * @return y
     * @throws IloException 
     */
    public IloNumVar IF_Between_Them_Y(String name, double lb, IloNumExpr exp, double ub) throws IloException {
        IloNumVar y = numVar(0, 1, name+".y");
        addIF_Between_Them_Y(name, y, lb, exp, ub);
        return y;
    }
    /**
     * <b>If</b> lb &le exp &le ub <b>then</b> y = 1;
     * @param name
     * @param y
     * @param lb
     * @param exp
     * @param ub
     * @throws IloException 
     */
    public void addIF_Between_Them_Y(String name, IloNumExpr y, double lb, IloNumExpr exp, double ub) throws IloException {
        IloIntVar z = boolVar("z");
        IloIntVar w = boolVar("w");
        addIF_Le_Them_Y(name, z, sum(exp, -ub));
        addIF_Ge_Them_Y(name, w, sum(exp, -lb));
        addAnd(name+".and", y, z, w);
    }
    
    private double print(double value){
        if(Double.isNaN(value)){
            return Double.NaN;
        }else if(Math.abs(value)<1e-12){
            return 0;
        }
        return value;
    }
    
    @Override
    public void print(String name, String format, double[] V) {
        //block(V);
        System.out.println("------------------------ [ "+name+"] ------------------------");
        //System.out.println("regulage_precision = "+regulage_presicion);
        for(double v : V){
            format(format, v);
            //System.out.printf(format+" ", print(v));
        }
        System.out.println();
        //unblock();
    }
    @Override
    public void print(String name, String format, double[][] V) {
        //block(V);
        System.out.println("------------------------ [ "+name+"] ------------------------");
        //System.out.println("regulage_precision = "+regulage_presicion);
        for (double[] v : V) {
            if(v!=null){
                for (int j = 0; j < v.length; j++) {
                    format(format, v[j]);
                    //System.out.printf(format+" ", print(v[j]));
                }
                System.out.println();
            }else{
                System.out.println("[null]");
            }
        }
        //unblock();
    }
    
    private void format(String format, double value){
        if(format.charAt(0) == '-'){
            String r = format.substring(format.indexOf(".")+1);
            int num = 0;
            while(r.length()>0 && r.charAt(0)<='9' && r.charAt(0)>='0'){
                num = 10*num + r.charAt(0)-'0';
                r = r.substring(1);
            }
            //System.out.println("format = ["+format+"] , num = " +num);
            
            if(Math.abs(value) < Math.pow(10, -num)){
                if(Math.abs(value)<Math.pow(10, -num-5)){
                    String s = String.format(format.substring(1), 0.0);
                    System.out.printf("%"+s.length()+"s ", "-");
                }else{
                    String s = String.format(format.substring(1), 0.0);
                    System.out.printf("%"+s.length()+"s ", "~");
                    //System.out.printf(format.substring(1)+" ", value);
                }
            }else{
                System.out.printf(format.substring(1)+" ", value);
                
//                if(Math.abs(value)<regulage_presicion){
//                    String s = String.format(format.substring(1), 0.0);
//                    System.out.printf("%"+s.length()+"s ", "-");
//                }else{
//                    System.out.printf(format.substring(1)+" ", value);
//                }
            }
        }else{
            System.out.printf(format+" ", value);
        }
    }
    
//    private boolean block = false;
//    private double regulage_presicion = 1e-9;
//    private void block(double[][]... V){
//        if(!block){
//            block = true;
//            double ub = Double.NEGATIVE_INFINITY;
//            double lb = Double.POSITIVE_INFINITY;
//            int count = 0;
//            for(double [][] v1:V){
//                for(double [] v2:v1){
//                    for(double v3:v2){
//                        ub = Math.max(ub, v3);
//                        lb = Math.min(lb, v3);
//                        count++;
//                    }
//                }
//            }
//            
//            regulage_presicion = lb + (ub-lb)/count ;
//            
//        }
//    }
//    private void block(double[] V){
//        if(!block){
//            block = true;
//            double ub = Double.NEGATIVE_INFINITY;
//            double lb = Double.POSITIVE_INFINITY;
//            int count = 0;
//            for(double v3:V){
//                ub = Math.max(ub, v3);
//                lb = Math.min(lb, v3);
//                count++;
//            }
//            
//            regulage_presicion = lb + (ub-lb)/count ;
//        }
//    }
//    private void unblock(){
//        if(block){
//            block = false;
//        }
//    }
    @Override
    public void print(String name, String format, double[][]... V) {
        //block(V);
        System.out.println("------------------------ [ "+name+"] ------------------------");
        //System.out.println("regulage_precision = "+regulage_presicion);
        int d2 = 0;
        for (double[][] v : V) {
            if(v!=null){
                d2 = Math.max(d2, v.length);
            }
        }
        
        for(int t=0; t<d2; t++){
            for (int n = 0; n < V.length; n++) {
                if(V[n] !=null && t<V[n].length && V[n][t]!=null){
                    if(t<V[n].length){
                        for (int j = 0; j < V[n][t].length; j++) {
                            format(format, V[n][t][j]);
                        }
                    }else{
                        for (int j = 0; j < V[n][0].length; j++) {
                            format(format, 0.0);
                        }
                    }
                    System.out.print("| ");
                }
            }
            System.out.println();
        }
        
        //unblock();
    }
    public void print(String name, String format, double[][] V, String formatCols, String ...cols) {
        //block(V);
        System.out.println("------------------------ [ "+name+"] ------------------------");
        //System.out.println("regulage_precision = "+regulage_presicion);
        for(String c : cols){
            System.out.printf(formatCols+" ", c);
        }
        System.out.println();
        for (double[] v : V) {
            if(v!=null){
                for (int j = 0; j < v.length; j++) {
                    format(format, v[j]);
                    //System.out.printf(format+" ", print(v[j]));
                }
                System.out.println();
            }else{
                System.out.println("[null]");
            }
        }
        //unblock();
    }
    /**
     * Deterministic, makes:
     * if Y them h*x^T le b
     * @param h
     * @param x
     * @param Y
     * @throws ilog.concert.IloException
     */
    public void addWaypointDet(Hyperplane h, IloNumVar x[], IloNumExpr Y) throws IloException {
        IloNumExpr exp = h.scalProd(this, x);
        exp = sum(exp, -h.b);
        addIF_Y_Them_Le("Ii", Y, exp);
    }
    /**
     * Deterministic, makes:
     * h*x^T le b true always
     * @param h
     * @param x
     * @throws ilog.concert.IloException
     */
    public void addWaypointDet(Hyperplane h, IloNumVar x[]) throws IloException {
        IloNumExpr exp = h.scalProd(this, x);
        exp = sum(exp, -h.b);
        addLe(exp, 0, "Ii");
    }

    /**
     * Deterministic, makes:
     * if Z them h*x^T ge b
     * @param h
     * @param x
     * @param Z
     * @throws ilog.concert.IloException
     */
    public void addObstacleDet(Hyperplane h, IloNumVar x[], IloNumExpr Z) throws IloException {
        IloNumExpr exp = h.scalProd(this, x);
        exp = sum(exp, -h.b);
        addIF_Y_Them_Ge("Oi", Z, exp);
    }
     /**
     * Uncertainty, makes:
     * if Y them h*x^T le b - c
     * @param h
     * @param x
     * @param Y
     * @param delta
     * @param uncertainty
     * @param DELTA
     * @param N
     * @throws ilog.concert.IloException
     * 
     */
    public void addWaypointUnc(Hyperplane h, IloNumVar x[], IloNumExpr Y, IloNumVar delta, double uncertainty, double DELTA, int N) throws IloException {
        if(uncertainty<1e-6){
            throw new IloException("uncertainty<1e-6");
        }else{
            //IloNumExpr risk = cplex.RiskAllocation(temp, uncertainty, approach.inst.chance_constraints[c], approach.Naprox(), "C(i,t,O)");
            IloNumExpr risk = RiskAllocation(delta, uncertainty, DELTA, N, "C(Oi,t)");
            
            IloNumExpr exp = h.scalProd(this, x);
            exp = sum(exp, -h.b);
            exp = sum(exp, prod(+1, risk));
            addIF_Y_Them_Le("Ii", Y, exp);
        }
    }
    /**
     * 
     * @param h
     * @param x
     * @param Y
     * @param erf_inv
     * @param uncertainty
     * @throws IloException 
     */
    public void addWaypointUnc(Hyperplane h, IloNumVar x[], IloNumExpr Y, IloNumExpr erf_inv, double uncertainty) throws IloException {
        if(uncertainty<1e-6){
            throw new IloException("uncertainty<1e-6");
        }else{
            //IloNumExpr risk = cplex.RiskAllocation(temp, uncertainty, approach.inst.chance_constraints[c], approach.Naprox(), "C(i,t,O)");
            IloNumExpr risk = RiskAllocation(erf_inv, uncertainty);
            
            IloNumExpr exp = h.scalProd(this, x);
            exp = sum(exp, -h.b);
            exp = sum(exp, prod(+1, risk));
            addIF_Y_Them_Le("Ii", Y, exp);
        }
    }
    /**
     * Uncertainty, makes:
     * h*x^T le b - c true always
     * @param h
     * @param x
     * @param delta
     * @param uncertainty
     * @param DELTA
     * @param N
     * @throws ilog.concert.IloException
     * 
     */
    public void addWaypointUnc(Hyperplane h, IloNumVar x[], IloNumVar delta, double uncertainty, double DELTA, int N) throws IloException {
        if(uncertainty<1e-6){
            throw new IloException("uncertainty<1e-6");
        }else{
            //IloNumExpr risk = cplex.RiskAllocation(temp, uncertainty, approach.inst.chance_constraints[c], approach.Naprox(), "C(i,t,O)");
            IloNumExpr risk = RiskAllocation(delta, uncertainty, DELTA, N, "C(Oi,t)");
            
            IloNumExpr exp = h.scalProd(this, x);
            exp = sum(exp, -h.b);
            exp = sum(exp, prod(+1, risk));
            addLe(exp, 0, "Ii");
        }
    }
    
    /**
     * Uncertainty, makes:
     * if Z them h*x^T ge b + c
     * @param h
     * @param x
     * @param Z
     * @param delta
     * @param uncertainty
     * @param DELTA
     * @param N
     * @throws ilog.concert.IloException
     * 
     */
    public void addObstacleUnc(Hyperplane h, IloNumVar x[], IloNumExpr Z, IloNumVar delta, double uncertainty, double DELTA, int N) throws IloException {
        if(uncertainty<1e-6){
            throw new IloException("uncertainty<1e-6");
        }else{
            //IloNumExpr risk = cplex.RiskAllocation(temp, uncertainty, approach.inst.chance_constraints[c], approach.Naprox(), "C(i,t,O)");
            IloNumExpr risk = RiskAllocation(delta, uncertainty, DELTA, N, "C(Oi,t)");
            
            IloNumExpr exp = h.scalProd(this, x);
            exp = sum(exp, -h.b);
            exp = sum(exp, prod(-1, risk));
            addIF_Y_Them_Ge("Oi", Z, exp);
        }
    }
    
    public void addObstacleUnc(Hyperplane h, IloNumVar x[], IloNumExpr Z, IloNumExpr erf_inv, double uncertainty) throws IloException {
        if(uncertainty<1e-6){
            throw new IloException("uncertainty<1e-6");
        }else{
            //IloNumExpr risk = cplex.RiskAllocation(temp, uncertainty, approach.inst.chance_constraints[c], approach.Naprox(), "C(i,t,O)");
            IloNumExpr risk = RiskAllocation(erf_inv, uncertainty);
            
            IloNumExpr exp = h.scalProd(this, x);
            exp = sum(exp, -h.b);
            exp = sum(exp, prod(-1, risk));
            addIF_Y_Them_Ge("Oi", Z, exp);
        }
    }

    /**
     * @param expr1
     * @param expr2
     * @return expr1 - expr2
     * @throws ilog.concert.IloException
     */
    public IloNumExpr minus(IloNumExpr expr1, IloNumExpr expr2) throws IloException {
       return sum(expr1, prod(-1,expr2));
    }
    public IloNumExpr minus(double value, IloNumExpr expr2) throws IloException {
       return sum(value, prod(-1,expr2));
    }
    public IloNumExpr minus(IloNumExpr expr1, double value) throws IloException {
       return sum(expr1, -value);
    }
    public IloNumExpr ScalProd(double[] a, IloNumVar[] x) throws IloException {
        IloNumExpr sum = null;
        int N = Math.min(a.length, x.length);
        for(int n=0; n<N; n++){
            sum = SumProd(sum, a[n], x[n]);
        }
        return sum;
    }
    /**
     * Calculate , makes:
     * if Z them h*x^T ge b + c
     * @param name
     * @param Zn
     * @param r
     * @param x
     * @param y 
     * @param MIN 
     * @param MAX
     * @param precision
     * @throws ilog.concert.IloException
     */
    public void aproxNorm2(String name, IloNumExpr[] Zn, IloNumExpr r, IloNumExpr x, IloNumExpr y, double MIN, double MAX, double precision) throws IloException {
        
//        //================================ calculating the BIG-M ============================
//        MAX += precision;
//        MIN -= Math.min(MIN,precision);
//        
//        double factor = 1.0/Math.cos(Math.PI/Zn.length) - 1; //factor of erro from poligonal of N faces to a perfect circle
//        double e_min = MIN*factor;                           //maximum error distance from the circle to the internal(MIN) poligon 
//        double e_max = MAX*factor;                           //maximum error distance from the circle to the external(MAX) poligon 
//        System.out.println("factor = "+factor);
//        System.out.println("e_min = "+e_min);
//        System.out.println("e_max = "+e_max);
//        double myBigM = MAX+MIN+e_max+e_min+3*precision+5;
//        //====================================================================================
            
        for(int n=0; n<Zn.length; n++){
            double cos = Math.cos(2*n*Math.PI/Zn.length);
            double sin = Math.sin(2*n*Math.PI/Zn.length);
//            IloNumExpr expr = sum(prod(cos, x), prod(sin, y), prod(-1, r) );
//            addGe(expr, prod(bigM, minus(Zn[n], 1)), name+".ge");
//            addLe(expr, 0, name+".le");
            
            //addIF_Y_Them_Ge(name+".ge", Zn[n], sum(prod(cos, x), prod(sin, y), prod(-1, r), constant(epsilon)));
            
            addGe(sum(prod(cos, x), prod(sin, y)), sum(-Math.min(MIN,precision), r, prod(bigM, minus(Zn[n], 1))), name+".ge");    //Ax - b >= M(y-1)  || -M + My
            addLe(sum(prod(cos, x), prod(sin, y)), sum(+precision, r), name+".le");
        }
    }
    
    public static void main(String[] args) throws IloException {
        CplexExtended cpx = new CplexExtended();
        IloNumVar Zn[] = cpx.boolVarArray(32, "Z");
        IloNumVar Vh = cpx.numVar(10, 30, "vh");
        IloNumVar Vx = cpx.numVar(Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, "vx");
        IloNumVar Vy = cpx.numVar(Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, "vy");
        
        IloNumVar Uh = cpx.numVar(1, 3, "uh");
        IloNumVar Ux = cpx.numVar(Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, "ux");
        IloNumVar Uy = cpx.numVar(Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, "uy");
        
        cpx.aproxNorm2("normV", Zn, Vh, Vx, Vy, 10, 30, 0.000001);
        cpx.aproxNorm2("normU", Zn, Uh, Ux, Uy, 1, 3, 0.000001);
        
        cpx.addEq(cpx.sum(Zn), 1);
        
        cpx.addMaximize(cpx.sum(cpx.minus(Vx, Ux),cpx.minus(Vy, Uy)));
        //cpx.addMaximize(cpx.sum(h,0));
        //cpx.addMinimize(cpx.sum(h,0));
        
        if(cpx.solve()){
            System.out.println("status = "+cpx.getStatus());
            System.out.println("objective = "+cpx.getObjValue());
            double alphaV = Math.atan2(cpx.getValue(Vy), cpx.getValue(Vx))*180/Math.PI;
            double alphaU = Math.atan2(cpx.getValue(Uy), cpx.getValue(Ux))*180/Math.PI;
            System.out.println("alphaV = "+alphaV);
            System.out.println("alphaU = "+alphaU);
            System.out.println("diff   = "+(alphaV-alphaU));
            
            System.out.println("V = ["+cpx.getValue(Vx)+" "+cpx.getValue(Vy)+"]");
            System.out.println("U = ["+cpx.getValue(Ux)+" "+cpx.getValue(Uy)+"]");
            System.out.println("V = ["+cpx.getValue(Vx)/cpx.getValue(Vh)+" "+cpx.getValue(Vy)/cpx.getValue(Vh)+"]");
            System.out.println("U = ["+cpx.getValue(Ux)/cpx.getValue(Uh)+" "+cpx.getValue(Uy)/cpx.getValue(Uh)+"]");
            
        }
    }

    public IloNumExpr sum(double value, IloNumExpr... expr) throws IloException {
        IloNumExpr sum = constant(value);
        for(IloNumExpr e : expr){
            sum = SumProd(sum, 1, e);
        }
        return sum;
    }

    /**<pre>
     * This encode the folowing lógic:
     *      z &ge &sum<sub>i</sub> |x<sub>i</sub>|
     * 
     * Thas is encoded as:
     *      z = &sum<sub>i</sub> (w<sub>i</sub>)    
     *      w<sub>i</sub> &ge +x<sub>i</sub>        &forall(i)
     *      w<sub>i</sub> &ge -x<sub>i</sub>        &forall(i)
     * Creates: 
     *      w<sub>i</sub> &ge 0
     * Complexity for i = 1 ... n
     *      B = 0       number of new binary variables
     *      C = n       number of new continous variables
     *      R = 2n+1    number of new restrictions
     * </pre>
     * @param name  is a simple name to identifie new variables or restrictions added by this call.
     * @param z     a expression for be defined as upper bound for the Manhattan norm. Any expression z &isin D &sube [-&infin +&infin] is valid.
     * @param Xi    a vector of expressions for Manhattan norm. Any expression x<sub>i</sub> &isin D<sub>i</sub>  &sube [-&infin +&infin] is valid.
     * @throws ilog.concert.IloException
     */
    public void addUpperBoundNormManhattan(String name, IloNumExpr z, IloNumExpr ...Xi) throws IloException {
        if(Xi.length<=1){
            addGe(z, prod(+1, Xi[0]), name+".pos");
            addGe(z, prod(-1, Xi[0]), name+".neg");
        }else{
            IloNumVar Wi[] = numVarArray(Xi.length, 0, Double.POSITIVE_INFINITY, name+".w");
            for(int i=0; i<Xi.length; i++){
                addGe(Wi[i], prod(+1, Xi[i]), name+".pos["+i+"]");
                addGe(Wi[i], prod(-1, Xi[i]), name+".neg["+i+"]");
            }
            addEq(z, sum(Wi), name+".sum");
        }
    }
    /**<pre>
     * This encode the folowing lógic:
     *      expr z;
     *      z &ge &sum<sub>i</sub> |x<sub>i</sub>|
     *      return z;
     * Thas is encoded as:
     *      var w<sub>i</sub> &ge 0;
     *      w<sub>i</sub> &ge +x<sub>i</sub>        &forall(i)
     *      w<sub>i</sub> &ge -x<sub>i</sub>        &forall(i)
     *      return &sum<sub>i</sub> (w<sub>i</sub>) 
     * Creates:
     *      w<sub>i</sub> &ge 0
     * Complexity for i = 1 ... n
     *      B = 0       number of new binary variables
     *      C = n       number of new continous variables
     *      R = 2n      number of new restrictions
     * </pre>
     * @param name  is a simple name to identifie new variables or restrictions added by this call.  
     * @param Xi    a vector of expressions for Manhattan norm. Any expression x<sub>i</sub> &isin D<sub>i</sub>  &sube [-&infin +&infin] is valid.
     * @return      a new expression z &ge 0, that is a upper bound for the Manhattan norm.
     * @throws ilog.concert.IloException
     */
    public IloNumExpr upperBoundNormManhattan(String name, IloNumExpr ...Xi) throws IloException {
        IloNumVar Wi[] = numVarArray(Xi.length, 0, Double.POSITIVE_INFINITY, name+".w");
        for(int i=0; i<Xi.length; i++){
            addGe(Wi[i], prod(+1, Xi[i]), name+".pos["+i+"]");
            addGe(Wi[i], prod(-1, Xi[i]), name+".neg["+i+"]");
        }
        return sum(Wi);
    }
//    
//    /**<pre>
//     * This encode the folowing aproximation:
//     *      z &#8819 x<sup>n</sup>                  for n &ge 1
//     * This aproximation is conservative once 
//     * it is ensure that:
//     *      z &ge x<sup>n</sup>
//     * It is encoded as:
//     *      z = &sum<sub>i</sub> (w<sub>i</sub>)    
//     *      w<sub>i</sub> &ge +x<sub>i</sub>        &forall(i)
//     *      w<sub>i</sub> &ge -x<sub>i</sub>        &forall(i)
//     * 
//     * Creates: 
//     *      w<sub>i</sub> &ge 0
//     * Complexity for i = 1 ... n
//     *      B = 0       number of new binary variables
//     *      C = n       number of new continous variables
//     *      R = 2n+1    number of new restrictions
//     * </pre>
//     * @param name
//     * @param z
//     * @param x
//     * @throws IloException 
//     */
//    public void addUpperBoundPower(String name, IloNumExpr z, IloNumExpr x, double n) throws IloException {
//        addUpperBoundPower(name, z, x, n);
//    }
//    /**<pre>
//     * 
//     * This encode the folowing aproximation:
//     *      z &#8819 f(x)                           &forall(x &isin D)
//     * This aproximation is conservative once 
//     * it is ensure that:
//     *      z &ge f(x)                              &forall(x &isin D)
//     * It is encoded as:
//     *      z = &sum<sub>i</sub> (w<sub>i</sub>)    
//     *      w<sub>i</sub> &ge +x<sub>i</sub>        &forall(i)
//     *      w<sub>i</sub> &ge -x<sub>i</sub>        &forall(i)
//     * 
//     * Creates: 
//     *      w<sub>i</sub> &ge 0
//     * Complexity for i = 1 ... n
//     *      B = 0       number of new binary variables
//     *      C = n       number of new continous variables
//     *      R = 2n+1    number of new restrictions
//     * </pre>
//     * @param name
//     * @param z
//     * @param x
//     * @throws IloException 
//     */
//    public void addUpperBoundConvexFunction(String name, IloNumExpr z, IloNumExpr x, iCplexConvexFunction f, double minX, double maxX) throws IloException {
//        addUpperBoundPower(name, z, x, n);
//    }
    
    
    /**
     * <b>If</b> constraints &ge 0 <b>then</b> y = 1;
     * @param name
     * @param y
     * @param constraints
     * @throws IloException 
     */
//    public final void addIF_Ge_Them_Y(String name, IloNumVar y, IloNumExpr ...constraints) throws IloException {
//        for(IloNumExpr exp : constraints){
//            addLe(exp, sum(-epsilon, prod(+bigM, y)), name);    //Ax - b <= My - e
//        }
//    }
    /**
     * y = 1 &lt==> constraints &ge 0
     * @param name
     * @param y
     * @param constraints
     * @throws IloException 
     */
//    public final void addEq_Y_Le(String name, IloNumVar y, IloNumExpr ...constraints) throws IloException {
//        addIF_Y_Them_Le(name, y, constraints);
//        addIF_Le_Them_Y(name, y, constraints);
//    }
    /**
     * y = 1 &lt==> constraints &ge 0
     * @param name
     * @param y
     * @param constraints
     * @throws IloException 
     */
//    public final void addEq_Y_Ge(String name, IloNumVar y, IloNumExpr ...constraints) throws IloException {
//        addIF_Y_Them_Ge(name, y, constraints);
//        addIF_Ge_Them_Y(name, y, constraints);
//    }
    
}
