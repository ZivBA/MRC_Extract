package meshi.energy.hydrogenBondsPlane;
import meshi.geometry.Distance;
import meshi.util.filters.Filter;

 public class CNtwoDistances{
     private Distance dis1 = null, dis2 = null;
     private byte counter = 0;
     private static final int TWO_DISTANCE_CAPACITY = 4;


     protected  CNtwoDistances() {}

     protected  CNtwoDistances (Distance dis1, Distance dis2) {
     this.dis1 = dis1;
     this.dis2 = dis2;
     counter = 2;
     }

     protected  int setDistance(Distance distance){
         if (counter == 0) {
             dis1 = distance;
             if (dis2 != null) throw new RuntimeException("Wrong order");
             counter++;
         }
         else if ((counter == 1) && (dis1 != distance)) {
             if (dis1 == null) throw new RuntimeException("Wrong order");
             if (dis2 != null) throw new RuntimeException("Couple dis2");             
             dis2 = distance;
             counter++;
         }
         return counter;
      }

     protected  void deleteDistance(Distance distance){
         if (counter == 0)
            throw new RuntimeException("try to delete from Empty Element");
         if (counter == 1)
         {
             dis1 = null;
             counter--;
             if (dis2 != null)
             throw new RuntimeException("empty distance1 but non empty distance2");
         }
         else if (counter == 2) {
             dis1 = dis2;
             dis2 = null;
             if (dis1 == null)
             throw new RuntimeException("empty distance1 but non empty distance2");
             counter--;
         }
      }

     protected Distance distance1(){return dis1;}
     protected Distance distance2(){return dis2;}
     protected Byte counter(){return counter;}
     protected boolean isInDM (){return (counter == 2);}

  /*
   public AtomList setAtoms(){
       atoms = new AtomList(TWO_DISTANCE_CAPACITY );
       atoms.fastAdd(dis1.atom1());
       atoms.fastAdd(dis1.atom2());
       atoms.fastAdd(dis2.atom1());
       atoms.fastAdd(dis2.atom2());
       return atoms;
   }
  */
     protected static class IsCNtwoDistances implements Filter {
     public boolean accept(Object obj) {
         return (obj instanceof CNtwoDistances);
     }
     }
 }
