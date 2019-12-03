file="";
num=0;       // must be initalized
geo="";       // must be initalized

//val=param1;      // param1 passed via -D on cmd-line
//echo(val,param1); // outputs 17,17


 module cutPuzzlePiece(a,filetouse){
 
     intersection(){
        name = str("fragmentsphere/fragment_", a, ".stl");
        scale(600)
       import(file=name, centre=true);
     
     
      //  rotate([-7,-30,10])
      //  translate([10,0,-110])

       // import("geometry/Frankvase.stl");
     
     
        //transformation for amber cup
        //rotate([-7,-30,10])
        //translate([10,0,-110])
         //transformation for Saltdeanpot
        //translate([2,-7,4])
        //scale(25)
        import(filetouse);
       //import("geometry/saltdeanpot_1K.stl");
        // sphere(100);
    }
 }
//echo(num);
//echo(file);
cutPuzzlePiece(num,file);

