/* This file is part of JGrasstools (http://www.jgrasstools.org)
 * (C) HydroloGIS - www.hydrologis.com 
 * 
 * JGrasstools is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package od;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import oms3.annotations.*;
import oms3.io.CSVTableWriter;
import org.jgrasstools.gears.libs.modules.JGTConstants;
import org.jgrasstools.gears.libs.modules.JGTModel;
import org.jgrasstools.gears.libs.monitor.IJGTProgressMonitor;
import org.jgrasstools.gears.libs.monitor.LogProgressMonitor;

@Description("Teorethical semivariogram models.")
@Documentation("vgm.html")
@Author(name = "Giuseppe Formetta, Adami Francesco", contact = " http://www.ing.unitn.it/dica/hp/?user=rigon")
@Keywords("Kriging, Hydrology")
@Label(JGTConstants.STATISTICS)
@Name("kriging")
@Status(Status.CERTIFIED)
@License("General Public License Version 3 (GPLv3)")
@SuppressWarnings("nls")
public class VGM extends JGTModel {
    
    @Description("Distances vector.")
    @In
    public double[] distnce = null;
    
    @Description("Sill value.")
    @In
    public double sill;
    
    @Description("Range value.")
    @In
    public double range;
    
    @Description("Nugget value.")
    @In
    public double nugget;
    
    @Description("Model name")
    @In
    public String modelname;
    
    @Description("The progress monitor.")
    @In
    public IJGTProgressMonitor pm = new LogProgressMonitor();
    
    @Description("Semivariance vector.")
    @Out
    public double[] result = null;
    
    @Description("Measured vector.")
    @Out
    public double[] obs = null;
    
    @Description("Measured vector.")
    @In
    public double[] inp = null;

    @Execute
    public void process() throws Exception {
        result = calculate(distnce, modelname, sill, range, nugget);
        obs = inp;
    }

//    public static double[][] vgm(double[] distance, String model, double Sill, double Range, double Nug) {
//        double[][] result = new double[distance.length][2];
//        double[] theor_variog = calculate(distance, model, Sill, Range, Nug);
//        for (int i = 0; i < distance.length; i++) {
//            result[i][0] = theor_variog[i];
//            result[i][1] = distance[i];
//            //System.out.println(result[i][0]);
//        }
//        return result;
//    }

    public static double[] calculate(double[] distance, String model, double sill, double range, double nug) {
        
        double[] result = null;

        if ("exponential".equals(model)) {
            result = fn_exponential(distance, sill, range, nug);
        }
        if ("gaussian".equals(model)) {
            result = fn_gaussian(distance, sill, range, nug);
        }
        if ("spherical".equals(model)) {
            result = fn_spherical(distance, sill, range, nug);
        }
        if ("pentaspherical".equals(model)) {
            result = fn_pentaspherical(distance, sill, range, nug);
        }
        if ("linear".equals(model)) {
            result = fn_linear(distance, sill, range, nug);
        }
        if ("circular".equals(model)) {
            result = fn_circolar(distance, sill, range, nug);
        }
        if ("bessel".equals(model)) {
            result = fn_bessel(distance, sill, range, nug);
        }
        if ("periodic".equals(model)) {
            result = fn_periodic(distance, sill, range, nug);
        }
        if ("hole".equals(model)) {
            result = fn_hole(distance, sill, range, nug);
        }
        if ("logaritmic".equals(model)) {
            result = fn_logaritmic(distance, sill, range, nug);
        }
        if ("power".equals(model)) {
            result = fn_power(distance, sill, range, nug);
        }
        if ("spline".equals(model)) {
            result = fn_spline(distance, sill, range, nug);
        }

//        System.out.println("Result + " + Arrays.toString(result));
        return result;
    }

    double[] fn_nugget(double[] dist, double range) {
        double[] nugget = new double[dist.length];
        for (int i = 0; i < dist.length; i++) {
            nugget[i] = (dist[i] == 0.0 ? 0.0 : 1.0);
        }
        return nugget;
    }

    public static double[] fn_exponential(double[] dist, double sill, double range, double nug) {
        int length = dist.length;
        double[] func = new double[length];
        for (int i = 0; i < length; i++) {
            if (dist[i] != 0.0) {
                func[i] = nug + sill * (1 - (Math.exp(-dist[i] / range)));
            }
            //System.out.println(func[i]);
        }
        return func;
    }

    public static double[] fn_gaussian(double[] dist, double sill, double range, double nug) {
        int length = dist.length;
        double[] func = new double[length];
        for (int i = 0; i < length; i++) {
            double hr;
            hr = dist[i] / (range);
            if (dist[i] != 0.0) {
                func[i] = nug + sill * (1.0 - (Math.exp(-(hr * hr))));
            }
            //System.out.println(func[i]);

        }
        return func;
    }

    public static double[] fn_spherical(double[] dist, double sill, double range, double nug) {
        int length = dist.length;
        double[] func = new double[length];
        for (int i = 0; i < length; i++) {
            double hr;
            hr = dist[i] / (range);
            if (dist[i] != 0.0) {
                func[i] = nug + sill * hr * (1.5 - 0.5 * hr * hr);
            }
            if (dist[i] >= range) {
                func[i] = sill;
            }
            //System.out.println(func[i]);
        }
        return func;

    }

    public static double[] fn_pentaspherical(double[] dist, double sill, double range, double nug) {
        int length = dist.length;
        double[] func = new double[length];
        double hr = 0, h2r2;
        for (int i = 0; i < length; i++) {
            hr = dist[i] / (range);
            h2r2 = hr * hr;
            if (dist[i] != 0.0) {
                func[i] = nug + sill * (hr * ((15.0 / 8.0) + h2r2 * ((-5.0 / 4.0) + h2r2 * (3.0 / 8.0))));
            }
            if (dist[i] >= range) {
                func[i] = sill;
            }
            //System.out.println(func[i]);
        }
        return func;

    }

    public static double[] fn_linear(double[] dist, double sill, double range, double nug) {
        int length = dist.length;
        double[] func = new double[length];
        for (int i = 0; i < length; i++) {
            if (dist[i] != 0.0) {
                func[i] = nug + sill * (dist[i] / range);
            }
            if (dist[i] >= range) {
                func[i] = sill;
            }
            //System.out.println(func[i]);
        }
        return func;

    }

    public static double[] fn_circolar(double[] dist, double sill, double range, double nug) {
        int length = dist.length;
        double hr;
        double[] func = new double[length];
        for (int i = 0; i < length; i++) {
            hr = dist[i] / (range);
            if (dist[i] != 0.0) {
                func[i] = nug + sill * ((2.0 / Math.PI) * (hr * Math.sqrt(1.0 - hr * hr) + Math.asin(hr)));
            }
            if (dist[i] >= range) {
                func[i] = sill;
            }
            //System.out.println(func[i]);

        }
        return func;

    }

    public static double[] fn_bessel(double[] dist, double sill, double range, double nug) {
        double hr;
        double MIN_BESS = 1.0e-3;
        double[] func = new double[dist.length];
        for (int i = 0; i < dist.length; i++) {
            hr = dist[i] / (range);
            if (hr > MIN_BESS) {
                func[i] = nug + sill * (1.0 - hr * bessk1(hr));
            }
            //System.out.println(func[i]);
        }
        return func;
    }

    static double bessk1(double x) /*
     * bessk1 from numerical recipes
     */ {
        double y, ans;

        if (x <= 2.0) {
            y = x * x / 4.0;
            ans = (Math.log(x / 2.0) * bessi1(x)) + (1.0 / x) * (1.0 + y * (0.15443144
                    + y * (-0.67278579 + y * (-0.18156897 + y * (-0.1919402e-1
                    + y * (-0.110404e-2 + y * (-0.4686e-4)))))));
        } else {
            y = 2.0 / x;
            ans = (Math.exp(-x) / Math.sqrt(x)) * (1.25331414 + y * (0.23498619
                    + y * (-0.3655620e-1 + y * (0.1504268e-1 + y * (-0.780353e-2
                    + y * (0.325614e-2 + y * (-0.68245e-3)))))));
        }
        return (double) ans;
    }

    static double bessi1(double x) /*
     * bessi1 from numerical recipes
     */ {
        double ax, ans;
        double y;

        if ((ax = Math.abs(x)) < 3.75) {
            y = x / 3.75;
            y *= y;
            ans = ax * (0.5 + y * (0.87890594 + y * (0.51498869 + y * (0.15084934
                    + y * (0.2658733e-1 + y * (0.301532e-2 + y * 0.32411e-3))))));
        } else {
            y = 3.75 / ax;
            ans = 0.2282967e-1 + y * (-0.2895312e-1 + y * (0.1787654e-1
                    - y * 0.420059e-2));
            ans = 0.39894228 + y * (-0.3988024e-1 + y * (-0.362018e-2
                    + y * (0.163801e-2 + y * (-0.1031555e-1 + y * ans))));
            ans *= (Math.exp(ax) / Math.sqrt(ax));
        }
        return (double) x < 0.0 ? -ans : ans;
    }

    public static double[] fn_periodic(double[] dist, double sill, double range, double nug) {
        double[] func = new double[dist.length];
        for (int i = 0; i < dist.length; i++) {
            if (dist[i] != 0.0) {
                func[i] = nug + sill * (1.0 - Math.cos(2.0 * Math.PI * dist[i] / (range)));
            }
            //System.out.println(func[i]);
        }
        return func;

    }

    public static double[] fn_hole(double[] dist, double sill, double range, double nug) {
        double[] func = new double[dist.length];
        for (int i = 0; i < dist.length; i++) {
            if (dist[i] != 0.0) {
                func[i] = nug + sill * (1.0 - Math.sin(dist[i] / (range)) / (dist[i] / (range)));
            }
            //System.out.println(func[i]);
        }
        return func;

    }

    public static double[] fn_logaritmic(double[] dist, double sill, double range, double nug) {
        double[] func = new double[dist.length];
        for (int i = 0; i < dist.length; i++) {
            if (dist[i] != 0.0) {
                func[i] = nug + sill * (Math.log(dist[i] / range));
            }
            //System.out.println(func[i]);
        }
        return func;

    }

    public static double[] fn_power(double[] dist, double sill, double range, double nug) {
        double[] func = new double[dist.length];
        for (int i = 0; i < dist.length; i++) {
            if (dist[i] != 0.0) {
                func[i] = nug + sill * (Math.pow(dist[i], range));
            }
            //System.out.println(func[i]);
        }
        return func;

    }

    public static double[] fn_spline(double[] dist, double sill, double range, double nug) {
        double[] func = new double[dist.length];
        for (int i = 0; i < dist.length; i++) {
            if (dist[i] != 0.0) {
                func[i] = nug + sill * (dist[i] * dist[i] * Math.log(dist[i]));
            }
            if (dist[i] >= range) {
                func[i] = sill;
            }
            //System.out.println(func[i]);
        }
        return func;

    }
    
    
    public static void main(String[] args) {
        String a = "gaussian";
        String b = "gaus" + "sian";
        
        System.out.println(a == b);
        System.out.println(a.equals(b));
    }
}
