FOLDER = "/joeschmoe/Desktop/"

def run_TDVP(L, XY, Z alpha):

    product_state = ["up"] * (L // 2) + ["down"] * (L - L // 2);
    psi = MPS.from_product_state(model1.lat.mps_sites(),
                                 product_state,
                                 bc=model1.lat.bc_MPS,
                                 form='B')

    LR_coeffs = long_range_coeffs(poly_decay, alpha, 1000, 10);
    Jx_ray = XY * LR_coeffs[1] * LR_coeffs[0];
    Jz_ray = Z * LR_coeffs[1] * LR_coeffs[0];
    xi_ray = LR_coeffs[0];

    model_params = dict(L=70, Jxx=Jx_ray, Jz=Jz_ray, xi=xi_ray, verbose=1)
    model1 = LongRangeHeisenberg(model_params)

    chi = 400;
    delta_t = 0.05
    tdvp_params = {
        'start_time': 0,
        'dt': delta_t,
        'trunc_params': {
            'chi_max': chi,
            'svd_min': 1.e-10,
            'trunc_cut': None
        }
    }

    tdvp_engine = tdvp.Engine(psi, model1, tdvp_params)
    times = [];
    Sent = [];
    Sz = [psi.expectation_value('Sz')];

    #Run 2-site TDVP to grow bond dimesion
    for i in range(200):
        print(i)
        tdvp_engine.run_two_sites(N_steps=1)
        times.append(tdvp_engine.evolved_time)
        Sent.append(psi.entanglement_entropy())
        Sz.append(psi.expectation_value('Sz'))

    #Run 1-site TDVP for most accurate long-time dynamics
    for i in range(500):
        print(i)
        tdvp_engine.run_one_site(N_steps=1)
        #psi_2=copy.deepcopy(psi)
        #psi_2.canonical_form()
        times.append(tdvp_engine.evolved_time)
        Sent.append(psi.entanglement_entropy())
        Sz.append(psi.expectation_value('Sz'))

    df_columns = ["alpha", "Entanglement Entropy", "Sz"];
    df = pd.DataFrame(columns=df_columns);
    df = df.append({"alpha": 0.5, "Entanglement Entropy": Sent, "Sz": Sz}, ignore_index=True);
    df.to_pickle(FOLDER + "/long_range_DW_alpha=%1.0dp%2.0dSz=%1.0dp%2.0d.pickle" % (int(alpha), int((alpha % 1)*100), int(Sz), int((Sz%1)*100)));
