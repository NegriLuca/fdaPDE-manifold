mesh_quadratounitario <- function(n){

    n_nodes = n*n;
    n_elements = (n-1)*(n-1);
    x = seq(0,1,length.out=n);
    y = seq(0,1,length.out=n);
    coord = matrix(0,n_nodes,2);
    k = 1;
    for(i in 1:n) {
        for(j in 1:n) {
            coord[k,] = c(x[i], y[j]);
            k = k + 1;
        }
    }
    seg = matrix(0,3*n_elements + 2*(n-1),2);
    k = 1;
    for(i in 1:(n-1)) {
        for(j in 1:(n-1)) {
            idx = j + (i-1)*n;
            seg[k,] = c(idx, idx+1);
            seg[k+1,] = c(idx, idx+n);
            seg[k+2,] = c(idx, idx+n+1);
            k = k+3;
        }
    }
    for(i in 1:(n-1)) {
        seg[k,] = c(n*i, n*(i+1));
        k = k+1;
    }
    for(i in 1:(n-1)) {
        #seg[k,] = c(n_nodes-n+1, n_nodes-n+i+1);
        seg[k,]= c(n_nodes-n+i,n_nodes-n+i+1);
        k = k+1;
    }

    meshf = list (nodes=coord, segments=seg)

    return(meshf)

}
