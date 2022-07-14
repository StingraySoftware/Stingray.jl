function positive_fft_bins(n_bin::Integer; include_zero::Bool = false)
    minbin = 2
    if include_zero
        minbin = 1
    end
    return (minbin : (n_bin+1) ÷ 2)
end

function poisson_level(norm::String; meanrate = nothing, n_ph = nothing, backrate::Real = 0.0)
    if norm == "abs"
        return 2.0 * meanrate
    elseif norm == "frac"
        return 2.0 / (meanrate - backrate)^2 * meanrate
    elseif norm == "leahy"
        return 2.0
    elseif norm == "none"
        return float(n_ph)
    else 
        throw(ArgumentError("Unknown value for norm: $norm"))
    end
end

function normalize_frac(unnorm_power::AbstractVector{<:Number}, dt::Real, n_bin::Integer, 
                        mean_flux::Real; background_flux::Real=0.0)
    if background_flux > 0
        power = @. unnorm_power * 2 * dt / ((mean_flux - background_flux) ^ 2 *
                                          n_bin)
    else
        # Note: this corresponds to eq. 3 in Uttley+14
        power = @. unnorm_power * 2 * dt / (mean_flux ^ 2 * n_bin)
    end
    return power
end

normalize_abs(unnorm_power::AbstractVector{<:Number}, dt::Real, n_bin::Integer) = 
    @. unnorm_power * 2 / n_bin / dt

normalize_leahy_from_variance(unnorm_power::AbstractVector{<:Number}, 
                              variance::Real, n_bin::Integer) = 
    @. unnorm_power * 2 / (variance * n_bin)

normalize_leahy_poisson(unnorm_power::AbstractVector{<:Number}, n_ph::Real) = 
    @. unnorm_power * 2 / n_ph

function normalize_periodograms(unnorm_power::AbstractVector{<:Number}, dt::Real, 
                                n_bin::Integer; mean_flux=nothing, n_ph=nothing,
                                variance=nothing, background_flux::Real=0.0, 
                                norm::String="frac",power_type::String="all")
    
    if norm == "leahy" && !isnothing(variance)
        pds = normalize_leahy_from_variance(unnorm_power, variance, n_bin)
    elseif norm == "leahy"
        pds = normalize_leahy_poisson(unnorm_power, n_ph)
    elseif norm == "frac"
        pds = normalize_frac(
            unnorm_power, dt, n_bin, mean_flux,
            background_flux=background_flux)
    elseif norm == "abs"
        pds = normalize_abs(unnorm_power, dt, n_bin)
    elseif norm == "none"
        pds = unnorm_power
    else 
        throw(ArgumentError("Unknown value for norm: $norm"))
    end

    if power_type == "all"
        return pds
    elseif power_type == "real"
        return real(pds)
    elseif power_type in ["abs", "absolute"]
        return abs.(pds)
    else 
        throw(ArgumentError("Unknown value for power_type: $power_type"))
    end

end

function bias_term(power1::Real, power2::Real, power1_noise::Real, 
                   power2_noise::Real, n_ave::Integer;
                   intrinsic_coherence::Real=1.0)
    
    if n_ave > 500
        return 0.0
    end
    return power1 * power2 - intrinsic_coherence * (power1 - power1_noise) * (power2 - power2_noise) / n_ave
end

function raw_coherence(cross_power::Number, power1::Real, power2::Real, 
                               power1_noise::Real, power2_noise::Real, 
                               n_ave::Integer; intrinsic_coherence::Real=1.0)

    bsq = bias_term(power1, power2, power1_noise, power2_noise, n_ave;
                    intrinsic_coherence=intrinsic_coherence)
    num = real(cross_power * conj(cross_power)) - bsq
    if num < 0
        num = real(cross_power * conj(cross_power))
    end
    den = power1 * power2
    return num / den
end

function estimate_intrinsic_coherence(cross_power::Complex, power1::Real,
                                              power2::Real, power1_noise::Real, 
                                              power2_noise::Real, n_ave::Integer)
    new_coherence = 1.0
    old_coherence = 0.0
    count = 0
    while (!(≈(new_coherence, old_coherence, atol=0.01)) && count< 40)
        old_coherence = new_coherence
        bsq = bias_term.(power1, power2, power1_noise, power2_noise,
                        n_ave, intrinsic_coherence=new_coherence)
        den = (power1 - power1_noise) * (power2 - power2_noise)
        num = real(cross_power * conj(cross_power)) - bsq
        if num < 0
            num = real(cross_power * conj(cross_power))
        end
        new_coherence = num / den
        count += 1
    end
    return new_coherence                                       
end

function error_on_averaged_cross_spectrum(cross_power:: AbstractVector{<:Complex}, 
                                          seg_power:: AbstractVector{<:Real}, 
                                          ref_power:: AbstractVector{<:Real}, n_ave::Integer,
                                          seg_power_noise::Real, ref_power_noise::Real;
                                          common_ref::Bool=false)

    if n_ave < 30
        @warn "n_ave is below 30. Please note that the error bars 
        on the quantities derived from the cross spectrum 
        are only reliable for a large number of averaged 
        powers."
    end
    
    two_n_ave = 2 * n_ave
    if common_ref
        Gsq = real.(cross_power .* conj(cross_power))
        bsq = bias_term.(seg_power, ref_power, seg_power_noise, ref_power_noise,
                        n_ave)
        frac = @. (Gsq - bsq) / (ref_power - ref_power_noise)
        power_over_2n = ref_power / two_n_ave

        # Eq. 18
        dRe = dIm = dG = @. NaNMath.sqrt(power_over_2n * (seg_power - frac))
        # Eq. 19
        dphi = @. NaNMath.sqrt(power_over_2n * (seg_power / (Gsq - bsq) -
                       1 / (ref_power - ref_power_noise)))

    else
        PrPs = ref_power .* seg_power
        dRe = @. NaNMath.sqrt((PrPs + real(cross_power) ^ 2 - imag(cross_power) ^ 2) /
                      two_n_ave)
        dIm = @. NaNMath.sqrt((PrPs - real(cross_power) ^ 2 + imag(cross_power) ^ 2) /
                      two_n_ave)
        gsq = raw_coherence.(cross_power, seg_power, ref_power,
                            seg_power_noise, ref_power_noise, n_ave)
        dphi = @. NaNMath.sqrt((1 - gsq) / (2 * gsq * n_ave))
        dG = sqrt.(PrPs ./ n_ave)
    end

    return dRe, dIm, dphi, dG
end

function cross_to_covariance(cross_power::Complex, ref_power::Real, 
                             ref_power_noise::Real, delta_nu::Real)
    # To be used with dot broadcast when need an array 
    return cross_power * sqrt(delta_nu / (ref_power - ref_power_noise))
end

function _which_segment_idx_fun(;binned::Bool=false, dt=nothing)
    if binned
        # Return a function, so that we can pass the correct dt as an argument.
        return (args...) -> generate_indices_of_segment_boundaries_binned(args...; dt=dt)
    else
        return generate_indices_of_segment_boundaries_unbinned 
    end
end

function get_average_ctrate(times:: AbstractVector{<:Real}, gti::AbstractMatrix{<:Real}, 
                            segment_size::Real; counts= nothing)
    n_ph = 0.0
    n_intvs = 0
    binned = !isnothing(counts)
    func = _which_segment_idx_fun(;binned)

    for (_, _, idx0, idx1) in func(times, gti, segment_size)
        if !(binned)
            n_ph += idx1 - idx0
        else
            n_ph += sum(@view counts[idx0+1:idx1])
        end
        n_intvs += 1
    end
    return (n_ph / (n_intvs * segment_size))
end

@resumable function get_flux_iterable_from_segments(times:: AbstractVector{<:Real}, 
                                                    gti::AbstractMatrix{<:Real}, 
                                                    segment_size::Real; n_bin=nothing,
                                                    fluxes=nothing, errors=nothing)
    dt = nothing
    binned = !isnothing(fluxes)
    if binned
        dt = Statistics.median(diff(@view times[1:100]))
    end
    fun = _which_segment_idx_fun(;binned, dt)

    for (s, e, idx0, idx1) in fun(times, gti, segment_size)
        if idx1 - idx0 < 2
            @yield nothing
            continue
        end
        if !binned
            event_times = @view times[idx0:idx1-1]
            cts = fit(Histogram,float.(event_times .- s);nbins=n_bin).weights
        else
            cts = float.(@view fluxes[idx0+1:idx1])
            if !isnothing(errors)
                cts = cts, @view errors[idx0+1:idx1]
            end
        end
        @yield cts
    end
end

function avg_pds_from_iterable(flux_iterable, dt::Real; norm::String="frac", 
                               use_common_mean::Bool=true,
                               silent::Bool=false)
    local_show_progress = show_progress
    if silent
        local_show_progress = identity
    end

    # Initialize stuff
    cross = unnorm_cross = nothing
    n_ave = 0
    freq = nothing

    sum_of_photons = 0
    common_variance = nothing
    fgt0 = nothing
    n_bin = nothing

    for flux in local_show_progress(flux_iterable)
        if isnothing(flux) || all(iszero,flux)
            continue
        end

        # If the iterable returns the uncertainty, use it to calculate the
        # variance.
        variance = nothing
        if flux isa Tuple
            flux, err = flux
            variance = Statistics.mean(err) ^ 2
        end

        # Calculate the FFT
        n_bin = length(flux)
        ft = fft(flux)

        # This will only be used by the Leahy normalization, so only if
        # the input light curve is in units of counts/bin
        n_ph = sum(flux)
        unnorm_power = real.(ft .* conj.(ft))

        # Accumulate the sum of means and variances, to get the final mean and
        # variance the end
        sum_of_photons += n_ph

        if !isnothing(variance)
            common_variance = sum_if_not_none_or_initialize(common_variance, variance)
        end

        # In the first loop, define the frequency and the freq. interval > 0
        if isnothing(cross)
            fgt0 = positive_fft_bins(n_bin)
            freq = fftfreq(n_bin, dt)[fgt0]
        end

        # No need for the negative frequencies
        keepat!(unnorm_power,fgt0)

        # If the user wants to normalize using the mean of the total
        # lightcurve, normalize it here
        cs_seg = unnorm_power
        if !(use_common_mean)
            mean = n_ph / n_bin

            cs_seg = normalize_periodograms(
                unnorm_power, dt, n_bin; mean_flux = mean, n_ph=n_ph,
                norm=norm, variance=variance,
            )
        end

        # Accumulate the total sum cross spectrum
        cross = sum_if_not_none_or_initialize(cross, cs_seg)
        unnorm_cross = sum_if_not_none_or_initialize(unnorm_cross,
                                                     unnorm_power)

        n_ave += 1
    end

    # If there were no good intervals, return nothing
    if isnothing(cross)
        return nothing
    end

    # Calculate the mean number of photons per chunk
    n_ph = sum_of_photons / n_ave
    # Calculate the mean number of photons per bin
    common_mean = n_ph / n_bin

    if !isnothing(common_variance)
        # Note: the variances we summed were means, not sums.
        # Hence M, not M*n_bin
        common_variance /= n_ave
    end

    # Transform a sum into the average
    unnorm_cross = unnorm_cross / n_ave
    cross = cross / n_ave

    # Final normalization (If not done already!)
    if use_common_mean
        cross = normalize_periodograms(
            unnorm_cross, dt, n_bin; mean_flux=common_mean, n_ph=n_ph,
            norm=norm, variance=common_variance
        )
    end

    results = DataFrame()
    results[!,"freq"] = freq
    results[!,"power"] = cross
    results[!,"unnorm_power"] = unnorm_cross
    results = attach_metadata(results,(n= n_bin, m= n_ave, dt= dt,
                         norm= norm,
                         df= 1 / (dt * n_bin),
                         nphots= n_ph,
                         mean= common_mean,
                         variance= common_variance,
                         segment_size= dt * n_bin))

    return results
end

function avg_cs_from_iterables_quick(flux_iterable1 ,flux_iterable2,
                                     dt::Real; norm::String="frac")
    unnorm_cross = unnorm_pds1 = unnorm_pds2 = nothing
    n_ave = 0
    fgt0 = n_bin = freq = nothing

    sum_of_photons1 = sum_of_photons2 = 0

    for (flux1, flux2) in zip(flux_iterable1, flux_iterable2)
        if isnothing(flux1) || isnothing(flux2) || all(iszero,flux1) || all(iszero,flux2)
            continue
        end

        n_bin = length(flux1)

        # Calculate the sum of each light curve, to calculate the mean
        n_ph1 = sum(flux1)
        n_ph2 = sum(flux2)

        # At the first loop, we define the frequency array and the range of
        # positive frequency bins (after the first loop, cross will not be
        # nothing anymore)
        if isnothing(unnorm_cross)
            freq = fftfreq(n_bin, dt)
            fgt0 = positive_fft_bins(n_bin)
        end

        # Calculate the FFTs
        ft1 = fft(flux1)
        ft2 = fft(flux2)

        # Calculate the unnormalized cross spectrum
        unnorm_power = ft1 .* conj.(ft2)

        # Accumulate the sum to calculate the total mean of the lc
        sum_of_photons1 += n_ph1
        sum_of_photons2 += n_ph2

        # Take only positive frequencies
        keepat!(unnorm_power,fgt0)

        # Initialize or accumulate final averaged spectrum
        unnorm_cross = sum_if_not_none_or_initialize(unnorm_cross,
                                                     unnorm_power)

        n_ave += 1
    end

    # If no valid intervals were found, return only `nothing`s
    if isnothing(unnorm_cross)
        return nothing
    end

    # Calculate the mean number of photons per chunk
    n_ph1 = sum_of_photons1 / n_ave
    n_ph2 = sum_of_photons2 / n_ave
    n_ph = sqrt(n_ph1 * n_ph2)
    # Calculate the mean number of photons per bin
    common_mean1 = n_ph1 / n_bin
    common_mean2 = n_ph2 / n_bin
    common_mean = n_ph / n_bin

    # Transform the sums into averages
    unnorm_cross ./= n_ave

    # Finally, normalize the cross spectrum (only if not already done on an
    # interval-to-interval basis)
    cross = normalize_periodograms(
        unnorm_cross,
        dt,
        n_bin;
        mean_flux = common_mean,
        n_ph=n_ph,
        norm=norm,
        variance=nothing,
        power_type="all",
    )

    # No negative frequencies
    freq = freq[fgt0]

    results = DataFrame()
    results[!,"freq"] = freq
    results[!,"power"] = cross
    results[!,"unnorm_power"] = unnorm_cross
    results = attach_metadata(results,(n= n_bin, m= n_ave, dt= dt,
                         norm= norm,
                         df= 1 / (dt * n_bin),
                         nphots= n_ph,
                         nphots1= n_ph1, nphots2= n_ph2,
                         variance= nothing,
                         mean= common_mean,
                         mean1= common_mean1,
                         mean2= common_mean2,
                         power_type= "all",
                         fullspec= false,
                         segment_size= dt * n_bin))

    return results
end

function avg_cs_from_iterables(
    flux_iterable1,
    flux_iterable2,
    dt::Real;
    norm::String="frac",
    use_common_mean::Bool=true,
    silent::Bool=false,
    fullspec::Bool=false,
    power_type::String="all",
    return_auxil::Bool=false)
    
    local_show_progress = show_progress
    if silent
        local_show_progress = (a) -> a
    end
    # Initialize stuff
    cross = unnorm_cross = unnorm_pds1 = unnorm_pds2 = pds1 = pds2 = nothing
    n_ave = 0
    fgt0 = n_bin = freq = nothing

    sum_of_photons1 = sum_of_photons2 = 0
    common_variance1 = common_variance2 = common_variance = nothing

    for (flux1, flux2) in local_show_progress(zip(flux_iterable1,
                                                flux_iterable2))
        if isnothing(flux1) || isnothing(flux2) || all(iszero,flux1) || all(iszero,flux2)
            continue
        end

        # Does the flux iterable return the uncertainty?
        # If so, define the variances
        variance1 = variance2 = nothing
        if flux1 isa Tuple
            flux1, err1 = flux1
            variance1 = Statistics.mean(err1) ^ 2
        end
        if flux2 isa Tuple
            flux2, err2 = flux2
            variance2 = Statistics.mean(err2) ^ 2
        end

        # Only use the variance if both flux iterables define it.
        if isnothing(variance1) || isnothing(variance2) 
            variance1 = variance2 = nothing
        else
            common_variance1 = sum_if_not_none_or_initialize(common_variance1,
                                                             variance1)
            common_variance2 = sum_if_not_none_or_initialize(common_variance2,
                                                             variance2)
        end

        n_bin = length(flux1)

        # At the first loop, we define the frequency array and the range of
        # positive frequency bins (after the first loop, cross will not be
        # nothing anymore)
        if isnothing(cross)
            freq = fftfreq(n_bin, dt)
            fgt0 = positive_fft_bins(n_bin)
        end

        # Calculate the FFTs
        ft1 = fft(flux1)
        ft2 = fft(flux2)

        # Calculate the sum of each light curve, to calculate the mean
        n_ph1 = sum(flux1)
        n_ph2 = sum(flux2)
        n_ph = sqrt(n_ph1 * n_ph2)

        # Calculate the unnormalized cross spectrum
        unnorm_power = ft1 .* conj.(ft2)

        # If requested, calculate the auxiliary PDSs
        if return_auxil
            unnorm_pd1 = real.(ft1 .* conj.(ft1))
            unnorm_pd2 = real.(ft2 .* conj.(ft2))
        end

        # Accumulate the sum to calculate the total mean of the lc
        sum_of_photons1 += n_ph1
        sum_of_photons2 += n_ph2

        # Take only positive frequencies unless the user wants the full
        # spectrum
        if !(fullspec)
            keepat!(unnorm_power,fgt0)
            if return_auxil
                keepat!(unnorm_pd1,fgt0)
                keepat!(unnorm_pd2,fgt0)
            end
        end

        cs_seg = unnorm_power
        if return_auxil
            p1_seg = unnorm_pd1
            p2_seg = unnorm_pd2
        end

        # If normalization has to be done interval by interval, do it here.
        if !(use_common_mean)
            mean1 = n_ph1 / n_bin
            mean2 = n_ph2 / n_bin
            mean = n_ph / n_bin
            variance = nothing

            if !isnothing(variance1)
                variance = sqrt(variance1 * variance2)
            end

            cs_seg = normalize_periodograms(
                unnorm_power, dt, n_bin; mean_flux = mean, n_ph=n_ph, norm=norm,
                power_type=power_type, variance=variance
            )

            if return_auxil
                p1_seg = normalize_periodograms(
                    unnorm_pd1, dt, n_bin; mean_flux = mean1, n_ph=n_ph1, norm=norm,
                    power_type=power_type, variance=variance1
                )
                p2_seg = normalize_periodograms(
                    unnorm_pd2, dt, n_bin; mean_flux = mean2, n_ph=n_ph2, norm=norm,
                    power_type=power_type, variance=variance2
            )
            end
        end
        # Initialize or accumulate final averaged spectra
        cross = sum_if_not_none_or_initialize(cross, cs_seg)
        unnorm_cross = sum_if_not_none_or_initialize(unnorm_cross,
                                                     unnorm_power)

        if return_auxil
            unnorm_pds1 = sum_if_not_none_or_initialize(unnorm_pds1,
                                                        unnorm_pd1)
            unnorm_pds2 = sum_if_not_none_or_initialize(unnorm_pds2,
                                                        unnorm_pd2)
            pds1 = sum_if_not_none_or_initialize(pds1, p1_seg)
            pds2 = sum_if_not_none_or_initialize(pds2, p2_seg)
        end

        n_ave += 1
    end

    # If no valid intervals were found, return only `nothing`s
    if isnothing(cross)
        return nothing
    end

    # Calculate the mean number of photons per chunk
    n_ph1 = sum_of_photons1 / n_ave
    n_ph2 = sum_of_photons2 / n_ave
    n_ph = sqrt(n_ph1 * n_ph2)

    # Calculate the common mean number of photons per bin
    common_mean1 = n_ph1 / n_bin
    common_mean2 = n_ph2 / n_bin
    common_mean = n_ph / n_bin

    if !isnothing(common_variance1)
        # Note: the variances we summed were means, not sums. Hence M, not M*N
        common_variance1 /= n_ave
        common_variance2 /= n_ave
        common_variance = sqrt(common_variance1 * common_variance2)
    end

    # Transform the sums into averages
    cross ./= n_ave
    unnorm_cross ./= n_ave
    if return_auxil
        unnorm_pds1 ./= n_ave
        unnorm_pds2 ./= n_ave
    end

    # Finally, normalize the cross spectrum (only if not already done on an
    # interval-to-interval basis)
    if use_common_mean
        cross = normalize_periodograms(
            unnorm_cross,
            dt,
            n_bin;
            mean_flux = common_mean,
            n_ph=n_ph,
            norm=norm,
            variance=common_variance,
            power_type=power_type,
        )
        if return_auxil
            pds1 = normalize_periodograms(
                unnorm_pds1,
                dt,
                n_bin;
                mean_flux = common_mean1,
                n_ph=n_ph1,
                norm=norm,
                variance=common_variance1,
                power_type=power_type,
            )
            pds2 = normalize_periodograms(
                unnorm_pds2,
                dt,
                n_bin;
                mean_flux = common_mean2,
                n_ph=n_ph2,
                norm=norm,
                variance=common_variance2,
                power_type=power_type,
            )
        end
    end
    # If the user does not want negative frequencies, don't give them
    if !(fullspec)
        freq = freq[fgt0]
    end

    results = DataFrame()
    results[!,"freq"] = freq
    results[!,"power"] = cross
    results[!,"unnorm_power"] = unnorm_cross

    if return_auxil
        results[!,"pds1"] = pds1
        results[!,"pds2"] = pds2
        results[!,"unnorm_pds1"] = unnorm_pds1
        results[!,"unnorm_pds2"] = unnorm_pds2
    end

    results = attach_metadata(results,(n= n_bin, m= n_ave, dt= dt,
                         norm= norm,
                         df= 1 / (dt * n_bin),
                         segment_size= dt * n_bin,
                         nphots= n_ph,
                         nphots1= n_ph1, nphots2= n_ph2,
                         countrate1= common_mean1 / dt,
                         countrate2= common_mean2 / dt,
                         mean= common_mean,
                         mean1= common_mean1,
                         mean2= common_mean2,
                         power_type= power_type,
                         fullspec= fullspec,
                         variance= common_variance,
                         variance1= common_variance1,
                         variance2= common_variance2))

    return results
    
end

function avg_pds_from_events(times:: AbstractVector{<:Real}, gti::AbstractMatrix{<:Real}, 
                             segment_size::Real, dt::Real; norm::String="frac",
                             use_common_mean::Bool=true, silent::Bool=false, 
                             fluxes=nothing, errors=nothing)
    if isnothing(segment_size)
        segment_size = max(gti) - min(gti)
    end
    n_bin = round(Int,segment_size / dt)
    dt = segment_size / n_bin

    flux_iterable = get_flux_iterable_from_segments(times, gti, segment_size;
                                                    n_bin, fluxes=fluxes,
                                                    errors=errors)
    cross = avg_pds_from_iterable(flux_iterable, dt, norm=norm,
                                  use_common_mean=use_common_mean,
                                  silent=silent)
    if !isnothing(cross)
        attach_metadata(cross,(gti=gti,))
    end
    return cross
    
end

function avg_cs_from_events(times1:: AbstractVector{<:Real}, times2:: AbstractVector{<:Real}, 
                            gti::AbstractMatrix{<:Real}, segment_size::Real, dt::Real; 
                            norm::String="frac", use_common_mean::Bool=true, 
                            fullspec::Bool=false, silent::Bool=false,
                            power_type::String="all", fluxes1=nothing, fluxes2=nothing,
                            errors1=nothing, errors2=nothing, return_auxil=false)
    if isnothing(segment_size) 
        segment_size = max(gti) - min(gti)
    end
    n_bin = round(Int, segment_size / dt)
    # adjust dt
    dt = segment_size / n_bin

    flux_iterable1 = get_flux_iterable_from_segments(
        times1, gti, segment_size; n_bin, fluxes=fluxes1, errors=errors1
    )
    flux_iterable2 = get_flux_iterable_from_segments(
        times2, gti, segment_size; n_bin, fluxes=fluxes2, errors=errors2
    )

    is_events = all(isnothing,(fluxes1, fluxes2, errors1,
                                                errors2))

    if (is_events
            && silent
            && use_common_mean
            && power_type == "all"
            && !fullspec
            && !return_auxil)
        results = avg_cs_from_iterables_quick(
            flux_iterable1,
            flux_iterable2,
            dt;
            norm=norm
        )

    else
        results = avg_cs_from_iterables(
            flux_iterable1,
            flux_iterable2,
            dt;
            norm=norm,
            use_common_mean=use_common_mean,
            silent=silent,
            fullspec=fullspec,
            power_type=power_type,
            return_auxil=return_auxil
        )
    end
    if !isnothing(results)
        attach_metadata(results,(gti=gti,))
    end
    return results
end
