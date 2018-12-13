 module cutPuzzlePiece(a){

 
 //for (a =[0:11]){

      //  a=11;
     
     intersection(){
        name = str("fragmentsphere/fragment_", a, ".stl");
        scale(100)

        import(file=name, centre=true);
     
     
        rotate([-7,-30,10])
        translate([10,0,-110])

        import("geometry/ambercuppoisson_100K.stl");
   }
    
 }
 

if (mode == 0) {
     cutPuzzlePiece(0);
}
if (mode == 1) {
     cutPuzzlePiece(1);
} else if (mode == 2) {
     cutPuzzlePiece(2);
}else if (mode == 3) {
     cutPuzzlePiece(3);
}else if (mode == 4) {
     cutPuzzlePiece(4);
}else if (mode == 5) {
     cutPuzzlePiece(5);
}else if (mode == 6) {
     cutPuzzlePiece(6);
}else if (mode == 7) {
     cutPuzzlePiece(7);
}else if (mode == 8) {
     cutPuzzlePiece(8);
}else if (mode == 9) {
     cutPuzzlePiece(9);
}else if (mode == 10) {
     cutPuzzlePiece(10);
}else if (mode == 11) {
     cutPuzzlePiece(11);

} 