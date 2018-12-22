
file="";
num=0;       // must be initalized
//val=param1;      // param1 passed via -D on cmd-line
//echo(val,param1); // outputs 17,17


 module cutPuzzlePiece(a,filetouse){
 
     intersection(){
        name = str("fragmentsphere/fragment_", a, ".stl");
        scale(600)
        rotate([-18,-135,-5])
        import(file=name, centre=true);
     
     
      //  rotate([-7,-30,10])
      //  translate([10,0,-110])

       // import("geometry/ambercuppoisson_100K.stl");
     
     
        translate([2,-7,4])
        scale(25)
        import(filetouse);
    
    }
 }
echo(num);
echo(file);
cutPuzzlePiece(num,file);