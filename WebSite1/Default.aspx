<%@ Page Language="c#" Debug="true" Explicit="True" %>

<%@ Import Namespace="System.Data" %>
<%@ Import Namespace="System.Data.OleDb" %>
<%@ Import Namespace="System.Data.OleDb" %>
<%@ Import Namespace="System.Drawing" %>
<%@ Import Namespace="System.Web.Mail" %>
<%@ Import Namespace="System.IO" %>
<%@ Import Namespace="Microsoft.CSharp" %>
<%@ Import Namespace="System.Drawing.Imaging" %>
<%@ Import Namespace="System.CodeDom.Compiler" %>
<%@ Import Namespace="System.Reflection" %>
<%@ Import Namespace="System.Collections" %>

<script language="C#" runat="server">
    
    public class tAlgebra
    {
        public tAlgebra()
        {
        }

        public static void Cross3(double ax, double ay, double az,
            double bx, double by, double bz,
            ref double outx, ref double outy, ref double outz)
        {
            outx = ay * bz - az * by;
            outy = az * bx - ax * bz;
            outz = ax * by - ay * bx;
        }

        public static void Reflect(double inx, double iny, double inz,
            double mirrorx, double mirrory, double mirrorz,
            ref double outx, ref double outy, ref double outz)
        {
            double perp1x = 0.0, perp1y = 0.0, perp1z = 0.0;
            double perp2x = 0.0, perp2y = 0.0, perp2z = 0.0;

            Cross3(inx, iny, inz, mirrorx, mirrory, mirrorz,
                ref perp1x, ref perp1y, ref perp1z);
            Normalize(ref perp1x, ref perp1y, ref perp1z);

            Cross3(perp1x, perp1y, perp1z, mirrorx, mirrory, mirrorz,
                ref perp2x, ref perp2y, ref perp2z);
            Normalize(ref perp2x, ref perp2y, ref perp2z);

            double a = mirrorx;
            double b = perp2x;
            double c = inx;

            double x = mirrory;
            double y = perp2y;
            double z = iny;

            double i = mirrorz;
            double j = perp2z;
            double k = inz;


            double n = 0.0, m = 0.0;

            double eps = 1.0E-5;

            if (Math.Abs(a) < eps)
            {
                if (Math.Abs(i) < eps)
                {
                    outx = -inx;
                    outy = iny;
                    outz = -inz;
                    return;
                }
                else
                {
                    double dn = (y - (x * j) / i);

                    if (Math.Abs(dn) < eps)
                    {
                        outx = -inx;
                        outy = iny;
                        outz = -inz;
                        return;
                    }

                    n = (z - (x * k) / i) / dn;
                    m = (k - (j * n)) / i;
                }
            }
            else
            {
                double dn = (y - (x * b) / a);

                if (Math.Abs(dn) < eps)
                {
                    outx = -inx;
                    outy = iny;
                    outz = -inz;
                    return;
                }

                n = (z - (x * c) / a) / dn;
                m = (c - (b * n)) / a;
            }

            double v1x = mirrorx;
            double v1y = mirrory;
            double v1z = mirrorz;

            double v2x = perp2x;
            double v2y = perp2y;
            double v2z = perp2z;

            v1x *= m;
            v1y *= m;
            v1z *= m;

            v2x *= n;
            v2y *= n;
            v2z *= n;

            outx = v1x - v2x;
            outy = v1y - v2y;
            outz = v1z - v2z;

            return;
        }

        public static double Dot3(double x1, double y1, double z1, double x2, double y2, double z2)
        {
            return ((x1 * x2) + (y1 * y2) + (z1 * z2));
        }

        public static double GetCosAngleV1V2(double v1x, double v1y, double v1z, double v2x, double v2y, double v2z)
        {
            // cos(t) = (v.w) / (|v|.|w|) = (v.w) / 1
            return Dot3(v1x, v1y, v1z, v2x, v2y, v2z);
        }

        public static double modv(double vx, double vy, double vz)
        {
            return System.Math.Sqrt(vx * vx + vy * vy + vz * vz);
        }

        public static bool Normalize(ref double vx, ref double vy, ref double vz)
        {
            double mod_v = tAlgebra.modv(vx, vy, vz);
            double eps = 1.0E-20;

            if (Math.Abs(mod_v) < eps)
                return true;

            vx = vx / mod_v;
            vy = vy / mod_v;
            vz = vz / mod_v;
            return false;
        }


    }
    public class tObject
    {
        public tObject()
        {
        }

        // Final illumination of a point (vertex) = ambient + diffuse + specular

        /*
        Ambient Light Contribution
        Ambient light = background light
        Light that is scattered by the environment
        Frequently assumed to be constant
        Very simple approximation of global illumination
        No direction: independent of light position, object
        orientation, observer�s position or orientation 
        */

        // ambient RGBA reflectance of the material default = (0.2, 0.2, 0.2, 1.0)
        public double ambientR, ambientG, ambientB, ambientA;

        /*
        Diffuse Light Calculation
        Need to decide how much light the object point receive
        from the light source � based on Lambert�s Law    
        Lambert�s law: the radiant energy D that a small surface
        patch receives from a light source is:
        D = I x cos (q)
        */
        // diffuse RGBA reflectance of the material default = (0.8, 0.8, 0.8, 1.0)
        public double diffuseR, diffuseG, diffuseB, diffuseA;

        // specular RGBA reflectance of the material default = (0.0, 0.0, 0.0, 1.0)
        // Fresnel factor, which depends on angle of incidence and refractive index 
        // of material (which in turn depends on wavelength of light) 
        // The overall effect is that the specular reflection does depend on the angle 
        // of incidence of the light (it is brighter for grazing angles of incidence), 
        // and it does reflect the material colour (more so if the light is nearly 
        // normal to surface)
        /*
        Specular light contribution
        The bright spot on the object
        The result of total reflection of
        the incident light in a concentrate 
        specular = Ks x I x cos(f)^n
         * */
        public double specularR, specularG, specularB, specularA;

        // RGBA emitted light intensity of the material default = (0.0, 0.0, 0.0, 1.0)
        public double emissionR, emissionG, emissionB, emissionA;

        // specifies the RGBA specular exponent of the material default = 0
        public double shininess;

        /*
        Illumination from each light:  (*sum all)
        Illum = ambient + diffuse + specular
        = Ka x I + Kd x I x ( cos q) + Ks x I x cos(f)^n
         q = angle between light vector and normal to surface
        f is the angle between reflection light and viewer direction
        n = shininess
        */

    }
    public class tSphere : tObject
    {
        public tSphere(double x, double y, double z, double r, double clr, double clg, double clb)
        {
            cx = x;
            cy = y;
            cz = z;
            radius = r;
            clR = clr;
            clG = clg;
            clB = clb;
        }
        public static double GetCoord(double i1, double i2, double w1, double w2, double p)
        {
            return ((p - i1) / (i2 - i1)) * (w2 - w1) + w1;
        }
        void Move(double vx, double vy, double vz)
        {
            cx += vx;
            cy += vy;
            cz += vz;
        }
        void MoveTo(double vx, double vy, double vz)
        {
            cx = vx;
            cy = vy;
            cz = vz;
        }
        void RotX(double angle)
        {
            double y = cy * System.Math.Cos(angle) - cz * System.Math.Sin(angle);
            double z = cy * System.Math.Sin(angle) + cz * System.Math.Cos(angle);
            cy = y;
            cz = z;
        }
        void RotY(double angle)
        {
            double x = cx * System.Math.Cos(angle) - cz * System.Math.Sin(angle);
            double z = cx * System.Math.Sin(angle) + cz * System.Math.Cos(angle);
            cx = x;
            cz = z;
        }
        public static double GetSphereIntersec(double cx, double cy, double cz, double radius,
                           double px, double py, double pz, double vx, double vy, double vz)
        {
            // x-xo 2 + y-yo 2 + z-zo 2 = r 2
            // x,y,z = p+tv 
            // At2 + Bt + C = 0


            double A = (vx * vx + vy * vy + vz * vz);
            double B = 2.0 * (px * vx + py * vy + pz * vz - vx * cx - vy * cy - vz * cz);
            double C = px * px - 2 * px * cx + cx * cx + py * py - 2 * py * cy + cy * cy + pz * pz - 2 * pz * cz + cz * cz - radius * radius;
            double D = B * B - 4 * A * C;
            double t = -1.0;
            if (D >= 0)
            {
                double t1 = (-B - System.Math.Sqrt(D)) / (2.0 * A);
                double t2 = (-B + System.Math.Sqrt(D)) / (2.0 * A);
                if (t1 < t2) t = t1; else t = t2;
            }

            return t;
        }
        public void getNormal(double x1, double y1, double z1, ref double vx1, ref double vy1, ref double vz1)
        {
            vx1 = x1 - cx;
            vy1 = y1 - cy;
            vz1 = z1 - cz;
        }
        public double cx, cy, cz, radius, clR, clG, clB;
    }

    private void Page_Load(object sender, System.EventArgs e)
    {
        Bitmap newBitmap = new Bitmap(200, 200, PixelFormat.Format32bppArgb);
        Graphics g = Graphics.FromImage(newBitmap);

        Pen blackPen = new Pen(Color.Black);
        Color clrBackground = Color.Black;
        g.FillRectangle(new SolidBrush(clrBackground), new Rectangle(0, 0, 200, 200));

        Rectangle rect = new Rectangle(0, 0, 200, 200);

        ///////////////////////////////////////
        System.Collections.ArrayList obj3dArrayList;
        obj3dArrayList = new System.Collections.ArrayList();
        tSphere sph1 = new tSphere(0.01, 0.001, 10, 200.0, 0.0, 0.0, 255.0);
        
        // ambient properties for the material   
        sph1.ambientR = 0.329412;
        sph1.ambientG = 0.223529;
        sph1.ambientB = 0.027451;

        // specular properties for the material   
        sph1.specularR = 0.992157;
        sph1.specularG = 0.941176;
        sph1.specularB = 0.807843;
        sph1.shininess = 27.8974;

        sph1.diffuseR = 0.780392;
        sph1.diffuseG = 0.568627;
        sph1.diffuseB = 0.113725;
        
        obj3dArrayList.Add(sph1);
        
        Graphics graphics = g;

        double px = 0.0;
        double py = 0.0;
        double pz = 600.0;

        double lpx = -500;// 200.0;
        double lpy = -500;// 200.0;
        double lpz = 400.0;


        double fMax = 320.0;

        for (int i = rect.Left; i <= rect.Right; i++)
        {
            double x = tSphere.GetCoord(rect.Left, rect.Right, -fMax, fMax, i);

            for (int j = rect.Top; j <= rect.Bottom; j++)
            {
                double y = tSphere.GetCoord(rect.Top, rect.Bottom, fMax, -fMax, j);

                double t = 1.0E10;

                double vx = x - px, vy = y - py, vz = -pz;

                double mod_v = tAlgebra.modv(vx, vy, vz);
                vx = vx / mod_v;
                vy = vy / mod_v;
                vz = vz / mod_v;

                bool bShadow = false;

                tSphere spherehit = null;

                for (int k = 0; k < (int)obj3dArrayList.Count; k++)
                {
                    tSphere sphn = (tSphere)obj3dArrayList[k];
                    double taux = tSphere.GetSphereIntersec(sphn.cx, sphn.cy, sphn.cz, sphn.radius, px, py, pz, vx, vy, vz);
                    if (taux < 0) continue;

                    if (taux > 0 && taux < t)
                    {
                        t = taux;
                        spherehit = sphn;
                    }
                }

                Color color = Color.FromArgb(0, 0, 0);

                if (spherehit != null)
                {
                    double intersx = px + t * vx,
                           intersy = py + t * vy,
                           intersz = pz + t * vz;

                    double normalX = intersx - spherehit.cx,
                           normalY = intersy - spherehit.cy,
                           normalZ = intersz - spherehit.cz;

                    double lvX = lpx - intersx,
                           lvY = lpy - intersy,
                           lvZ = lpz - intersz;

                    tAlgebra.Normalize(ref normalX, ref normalY, ref normalZ);
                    tAlgebra.Normalize(ref lvX, ref lvY, ref lvZ);

                    double cost = tAlgebra.GetCosAngleV1V2(lvX, lvY, lvZ,
                        normalX, normalY, normalZ);

                    double cosf = 0;

                    double vReflX = 0, vReflY = 0, vReflZ = 0;

                    tAlgebra.Reflect(-lvX, -lvY, -lvZ,
                                     normalX, normalY, normalZ,
                                     ref vReflX, ref vReflY, ref vReflZ);

                    tAlgebra.Normalize(ref vReflX, ref vReflY, ref vReflZ);
                    tAlgebra.Normalize(ref vx, ref vy, ref vz);

                    cosf = tAlgebra.GetCosAngleV1V2(vx, vy, vz, vReflX, vReflY, vReflZ);

                    double result1 = Math.Max(0, cost) * 255.0;
                    double result2 = Math.Pow(Math.Max(0, cosf), spherehit.shininess) * 255.0;

                    double rgbR = (spherehit.ambientR * 255.0) + (spherehit.diffuseR * result1) + (spherehit.specularR * result2);
                    double rgbG = (spherehit.ambientG * 255.0) + (spherehit.diffuseG * result1) + (spherehit.specularG * result2);
                    double rgbB = (spherehit.ambientB * 255.0) + (spherehit.diffuseB * result1) + (spherehit.specularB * result2);

                    rgbR = Math.Min(rgbR, 255);
                    rgbG = Math.Min(rgbG, 255);
                    rgbB = Math.Min(rgbB, 255);
                    rgbR = Math.Max(0, rgbR);
                    rgbG = Math.Max(0, rgbG);
                    rgbB = Math.Max(0, rgbB);

                    color = Color.FromArgb((int)rgbR, (int)rgbG, (int)rgbB);
                }

                Brush brs = new SolidBrush(color);
                graphics.FillRectangle(brs, i, j, 1, 1);
                brs.Dispose();

            }// for pixels lines
        }// for pixels columns
        ///////////////////////////////////////

        MemoryStream tempStream = new MemoryStream();
        newBitmap.Save(tempStream, ImageFormat.Png);
        Response.ClearContent();
        Response.ContentType = "image/png";
        Response.BinaryWrite(tempStream.ToArray());
        Response.Flush();

    }


</script>

<html xmlns="http://www.w3.org/1999/xhtml">
<head>
    <title>Raytracing</title>
</head>
<body style="padding-right: 0px; padding-left: 0px; padding-bottom: 0px; margin: 0px;
    padding-top: 0px; text-align: center;">
    <form id="form1" runat="server">
    </form>

    <script src="http://www.google-analytics.com/urchin.js" type="text/javascript">
    </script>

    <script type="text/javascript">
_uacct = "UA-2275356-1";
urchinTracker();
    </script>

</body>
</html>
