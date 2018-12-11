//intersection(){
/*translate([-50,-70,100])
scale(100)

%import("fragmentsphere/fragment_0.stl");
%import("fragmentsphere/fragment_1.stl");
%import("fragmentsphere/fragment_2.stl");
%import("fragmentsphere/fragment_3.stl");
%import("fragmentsphere/fragment_4.stl");
%import("fragmentsphere/fragment_5.stl");
%import("fragmentsphere/fragment_6.stl");
rotate([-7,-30,10])
        translate([10,0,-100])
        
        import("geometry/ambercuppoisson_100K.stl");*/

//}

 for (a =[0:2]){


        name = str("fragmentsphere/fragment_", a, ".stl");
        scale(100)

        %import(file=name, centre=true);
     
     
        rotate([-7,-30,10])
        translate([10,0,-110])

        import("geometry/ambercuppoisson_100K.stl");

    
    
    
    }
