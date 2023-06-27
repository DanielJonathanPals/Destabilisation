using Random
using Kronecker

@testset "VAR model" begin

    B_hat = [-7.87301e-10  -4.91341e-8  1.0  8.67039e-10  2.82692e-8  1.0]
    Σ_hat_u = reshape([1.0728884458060339e-20],1,1)
    Σ_tilde_u = reshape([1.0664511151311976e-20],1,1)
    Σ_β_hat = [6.5246e-19    3.93324e-18  -8.21302e-22  -7.24205e-19   8.73169e-19  -4.34374e-18
                3.93324e-18   6.87154e-16   2.52076e-19  -4.3523e-18    7.42029e-17  -7.62544e-16
                -8.21302e-22   2.52076e-19   8.53403e-19   8.51485e-22   8.73971e-19  -3.67077e-19
                -7.24205e-19  -4.3523e-18    8.51485e-22   8.03858e-19  -1.01829e-18   4.80643e-18
                8.73169e-19   7.42029e-17   8.73971e-19  -1.01829e-18   5.41495e-16  -8.23982e-17
                -4.34374e-18  -7.62544e-16  -3.67077e-19   4.80643e-18  -8.23982e-17   8.46229e-16]
    Γ_hat = [ 1.0          0.00290959    0.000254886  0.900952      8.5257e-5    0.00264608
            0.00290959   0.000925389   8.5257e-5    0.00264608    2.21097e-7   0.000833841
            0.000254886  8.5257e-5     2.04489e-5   0.000231854  -6.06721e-10  7.68261e-5
            0.900952     0.00264608    0.000231854  0.811729      7.68261e-5   0.00240615
            8.5257e-5    2.21097e-7   -6.06721e-10  7.68261e-5    2.7427e-8    2.03174e-7
            0.00264608   0.000833841   7.68261e-5   0.00240615    2.03174e-7   0.000751363]
    Σ_x1_hat = reshape([1.07932577648087e-20],1,1)
    h(x,p) = x[1] .* p

    model = VARmodel(h,1,1,2,2,5,B_hat,Σ_hat_u,Σ_tilde_u,Σ_β_hat,Γ_hat,Σ_x1_hat)
    x_tr = reshape([-0.0113248],1,1)
    p_tr = reshape([-0.00656866,  0.908299],2,1)
    
    h1(x,p) = reshape(Array(x ⊗ p),length(x ⊗ p))
    h2(x) = x.^2

    model2 = VARmodel(h1,4,2,3,6,11,ones(2,45),ones(2,2),ones(2,2),ones(90,90),ones(45,45),ones(2,2))
    
    @test oneStepPred(model,x_tr,p_traj=p_tr) ≈ -0.016854963956433575 .* ones(1,1) atol=1e-4

    @test_throws ErrorException oneStepPred(model,ones(2,3),p_traj=ones(3,3))
    @test_throws ErrorException oneStepPred(model2,ones(2,3),p_traj=ones(3,3))



    @test_throws ErrorException VARmodel(nothing,4,2,3,6,11,ones(2,45),ones(2,2),ones(2,2),ones(90,90),ones(45,45),ones(2,2)) 
    @test_throws ErrorException VARmodel(h2,4,2,0,1,3,ones(2,13),ones(2,2),ones(2,2),ones(26,26),ones(13,13),ones(2,2)) 
    @test_throws ErrorException VARmodel(h1,4,2,3,5,10,ones(2,41),ones(2,2),ones(2,2),ones(82,82),ones(41,41),ones(2,2)) 
    @test_throws ErrorException VARmodel(h1,0,2,3,6,11,ones(2,45),ones(2,2),ones(2,2),ones(90,90),ones(45,45),ones(2,2))
    @test_throws ErrorException VARmodel(h1,4,2,3,6,10,ones(2,45),ones(2,2),ones(2,2),ones(90,90),ones(45,45),ones(2,2))
    @test_throws ErrorException VARmodel(h1,4,2,3,6,11,ones(2,44),ones(2,2),ones(2,2),ones(90,90),ones(45,45),ones(2,2))
    @test_throws ErrorException VARmodel(h1,4,2,3,6,11,ones(2,45),ones(2,3),ones(2,2),ones(90,90),ones(45,45),ones(2,2))
    @test_throws ErrorException VARmodel(h1,4,2,3,6,11,ones(2,45),ones(2,2),ones(1,2),ones(90,90),ones(45,45),ones(2,2))
    @test_throws ErrorException VARmodel(h1,4,2,3,6,11,ones(2,45),ones(2,2),ones(2,2),ones(91,90),ones(45,45),ones(2,2))
    @test_throws ErrorException VARmodel(h1,4,2,3,6,11,ones(2,45),ones(2,2),ones(2,2),ones(90,90),ones(44,45),ones(2,2))
    @test_throws ErrorException VARmodel(h1,4,2,3,6,11,ones(2,45),ones(2,2),ones(2,2),ones(90,90),ones(45,45),ones(1,2))
end