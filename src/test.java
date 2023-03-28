import jdk.incubator.vector.*;
public class test {
    public static void main(String[] args){
        float[] x = {2,3,4,5};
        float[] y = {3,4,5,6};
        float[] cCol = {7,8,9,10};
        var cy = FloatVector.fromArray(FloatVector.SPECIES_PREFERRED, cCol, 0);
       var vx = FloatVector.fromArray(FloatVector.SPECIES_PREFERRED, x, 0);
        var vy = FloatVector.fromArray(FloatVector.SPECIES_PREFERRED, y, 0);
        var vxsq = vx.mul(vx);
        var vysq = vy.mul(vy);
        var z = vxsq.add(vysq);
        var zy = cy.add(vx.mul(vy.mul(2)));
        float[] zs = new float[4];
        float[] zCol = new float[4];
        zy.intoArray(zCol, 0);
        z.intoArray(zs, 0);
        for ( int i = 0; i < zCol.length; i++){
            System.out.println(zCol[i] + " ");
        }
    }
}
