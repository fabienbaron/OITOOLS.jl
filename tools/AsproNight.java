// AsproNight — ASPRO 2's OBSERVABILITY for one night, printed as parseable lines.
//
// The companion to `AsproNoise.java`, which compares error bars; this one compares the
// night itself: when a target is observable, when it transits, and where the Sun puts the
// twilight boundaries. Same boot requirements — see AsproNoise.java's header for why the
// app has to be launched and why a display must exist.
//
//   javac -cp Aspro2-26.09.jar AsproNight.java
//   xvfb-run -a java -cp .:Aspro2-26.09.jar AsproNight --target Vega \
//       --ra 18:36:56.336 --dec +38:47:01.28 --date 2027-04-24 --minElev 30
//
// `--pops 451215` holds the PoPs (one digit per station, in the order of `--stations`) instead
// of letting ASPRO choose. `--date` is the night ENDING on that date's morning: the night that
// `night_plan` calls 2026-09-15 is ASPRO's 2026-09-16.
//
// Output is one `KEY<TAB>...` line per fact, so the Julia side can read it without parsing
// prose. Times are UTC, ISO 8601, because that is the one clock both codes agree on without
// a time-zone database.
//
// PINNED AGAINST Aspro2 26.09, same as AsproNoise.java. The signatures used here:
//   ObservabilityService(ObservationSetting).compute()
//   ObservabilityData: getMapStarVisibilities / getSunIntervals / getMoonIllumPercent
//   StarObservabilityData: getVisible / getTransitDate / getType / getTargetName
//   DateTimeInterval: getStartDate / getEndDate
//   SunTimeInterval: getStartDate / getEndDate / getType

import fr.jmmc.jmcs.Bootstrapper;
import fr.jmmc.aspro.Aspro2;
import fr.jmmc.aspro.model.ConfigurationManager;
import fr.jmmc.aspro.model.ObservationManager;
import fr.jmmc.aspro.model.oi.ObservationSetting;
import fr.jmmc.aspro.model.oi.Target;
import fr.jmmc.aspro.model.observability.ObservabilityData;
import fr.jmmc.aspro.model.observability.StarObservabilityData;
import fr.jmmc.aspro.model.observability.DateTimeInterval;
import fr.jmmc.aspro.model.observability.SunTimeInterval;
import fr.jmmc.aspro.service.ObservabilityService;
import fr.jmmc.aspro.model.TimeRef;
import fr.jmmc.aspro.model.observability.SunTimeInterval.SunType;
import fr.jmmc.aspro.Preferences;

import java.text.SimpleDateFormat;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TimeZone;

public class AsproNight {

    static Map<String, String> ARGS = new HashMap<String, String>();
    static SimpleDateFormat ISO = new SimpleDateFormat("yyyy-MM-dd'T'HH:mm:ss");

    static String s(String k, String dflt) {
        String v = ARGS.get(k);
        return (v == null) ? dflt : v;
    }

    static double d(String k, double dflt) {
        String v = ARGS.get(k);
        return (v == null) ? dflt : Double.parseDouble(v);
    }

    public static void main(String[] argv) throws Exception {
        for (int i = 0; i < argv.length; i++) {
            if (argv[i].startsWith("--")) {
                String key = argv[i].substring(2);
                String val = (i + 1 < argv.length && !argv[i + 1].startsWith("--")) ? argv[++i] : "true";
                ARGS.put(key, val);
            }
        }
        ISO.setTimeZone(TimeZone.getTimeZone("UTC"));
        Bootstrapper.launchApp(new Aspro2(new String[]{}));

        // ASPRO logs the per-beam optical path and the per-baseline W range at DEBUG. Those are
        // the only way to compare a delay model term by term rather than by its output.
        if (ARGS.containsKey("debug")) {
            ch.qos.logback.classic.Logger lg = (ch.qos.logback.classic.Logger)
                    org.slf4j.LoggerFactory.getLogger("fr.jmmc.aspro.service.ObservabilityService");
            lg.setLevel(ch.qos.logback.classic.Level.DEBUG);
        }

        ConfigurationManager cm = ConfigurationManager.getInstance();
        String interf   = s("interferometer", "CHARA");
        String conf     = s("config", cm.getInterferometerConfigurationNames(interf).lastElement());
        String instrum  = s("instrument", "MIRCX-MYSTIC");
        String stations = s("stations", "S1 S2 E1 E2 W1 W2");
        String mode     = s("mode", "Low_H");
        String name     = s("target", "Vega");

        ObservationManager om = ObservationManager.getInstance();
        om.reset();
        om.setInterferometerConfigurationName(conf);
        om.setInstrumentConfigurationName(instrum);
        om.setInstrumentConfigurationStations(new Object[]{stations});
        om.setInstrumentMode(mode);
        om.setWhen(new SimpleDateFormat("yyyy-MM-dd").parse(s("date", "2027-04-24")));
        om.setInstrumentSamplingPeriod(Double.valueOf(d("sampling", 40.0)));
        om.setInstrumentAcquisitionTime(Double.valueOf(d("acqTime", 600.0)));
        om.setMinElevation(d("minElev", 30.0));
        om.setAtmosphereQuality(fr.jmmc.aspro.model.oi.AtmosphereQuality.fromValue(s("atm", "Average")));
        // PoPs held instead of searched: one digit per station, in the order of --stations.
        if (ARGS.containsKey("pops")) {
            om.setInstrumentConfigurationPoPs(s("pops", ""));
        }

        Target t = new Target();
        t.setName(name);
        t.setRA(s("ra", "18:36:56.336"));
        t.setDEC(s("dec", "+38:47:01.28"));
        t.setFLUXV(Double.valueOf(d("magV", 2.0)));
        t.setFLUXH(Double.valueOf(d("magH", 2.0)));
        t.setFLUXK(Double.valueOf(d("magK", 2.0)));

        ObservationSetting obs = om.getMainObservation();
        obs.getTargets().add(t);

        // `new ObservabilityService(obs)` defaults to TimeRef.LST -- every Date it returns is
        // then a SIDEREAL time wearing a civil-time type, which reads as a seven-hour error if
        // you assume otherwise. The full constructor is the only way to ask for UTC.
        TimeRef ref = TimeRef.valueOf(s("timeref", "UTC").toUpperCase());
        ObservabilityData od = new ObservabilityService(obs, ref, false, false, false,
                SunType.Night,
                Preferences.getInstance().getBestPopsAlgorithm(),
                Preferences.getInstance().getBestPopsCriteriaSigma(),
                Preferences.getInstance().getBestPopsCriteriaAverageWeight()
        ).compute();
        System.out.println("TIMEREF\t" + od.getTimeRef());

        System.out.println("CONFIG\t" + conf + "\t" + instrum + "\t" + stations
                           + "\tminElev=" + d("minElev", 30.0));
        System.out.println("MOONILLUM\t" + od.getMoonIllumPercent());
        if (od.getDateMidnight() != null) {
            System.out.println("MIDNIGHT\t" + ISO.format(od.getDateMidnight()));
        }
        List<SunTimeInterval> suns = od.getSunIntervals();
        if (suns != null) {
            for (SunTimeInterval si : suns) {
                System.out.println("SUN\t" + si.getType() + "\t" + ISO.format(si.getStartDate())
                                   + "\t" + ISO.format(si.getEndDate()));
            }
        }
        // Which POP combination ASPRO settled on. Without it a delay-limited window cannot be
        // compared: the same target on the same night gives a different window per combination,
        // and the station ORDER in the configuration string is what a POP index refers to.
        if (od.getUsedPops() != null) {
            System.out.println("USEDPOPS\t" + od.getUsedPops().getIdentifier());
        }
        StringBuilder sb = new StringBuilder();
        if (od.getStationNameWithPops(sb)) {
            System.out.println("STATIONPOPS\t" + sb.toString());
        }
        if (od.getBestPopList() != null) {
            int k = 0;
            for (fr.jmmc.aspro.model.observability.PopCombination pc : od.getBestPopList()) {
                System.out.println("BESTPOP\t" + (++k) + "\t" + pc.getIdentifier());
                if (k >= 5) break;
            }
        }
        Map<String, List<StarObservabilityData>> vis = od.getMapStarVisibilities();
        if (vis != null) {
            for (Map.Entry<String, List<StarObservabilityData>> e : vis.entrySet()) {
                for (StarObservabilityData so : e.getValue()) {
                    if (so.getTransitDate() != null) {
                        System.out.println("TRANSIT\t" + e.getKey() + "\t" + so.getType()
                                           + "\t" + ISO.format(so.getTransitDate()));
                    }
                    // Azimuth and elevation along the track, so our own alt/az can be checked
                    // against ASPRO's rather than only the windows they produce.
                    if (so.getTargetPositions() != null) {
                        for (Map.Entry<java.util.Date,
                                fr.jmmc.aspro.model.observability.TargetPositionDate> tp
                                : so.getTargetPositions().entrySet()) {
                            System.out.println("POS\t" + e.getKey() + "\t"
                                    + ISO.format(tp.getKey()) + "\t"
                                    + tp.getValue().getHa() + "\t"
                                    + tp.getValue().getAzimuth() + "\t"
                                    + tp.getValue().getElevation());
                        }
                    }
                    // ASPRO draws TWO windows: `getVisible` honours the SOFT limits and is the
                    // solid bar; `getVisibleNoSoftLimits` is the dashed extension beyond them.
                    // For CHARA the soft limit is the 80 deg "motors report a warning" elevation,
                    // against the 85 deg hard maximum.
                    if (so.getVisibleNoSoftLimits() != null) {
                        for (DateTimeInterval in : so.getVisibleNoSoftLimits()) {
                            System.out.println("SOFTEXT\t" + e.getKey() + "\t"
                                    + ISO.format(in.getStartDate()) + "\t"
                                    + ISO.format(in.getEndDate()));
                        }
                    }
                    if (so.getVisibleNoSoftLimitsDesc() != null) {
                        System.out.println("SOFTDESC\t" + so.getVisibleNoSoftLimitsDesc());
                    }
                    List<DateTimeInterval> iv = so.getVisible();
                    if (iv == null) continue;
                    for (DateTimeInterval in : iv) {
                        System.out.println("VISIBLE\t" + e.getKey() + "\t" + so.getType()
                                           + "\t" + so.getLegendLabel(0)
                                           + "\t" + ISO.format(in.getStartDate())
                                           + "\t" + ISO.format(in.getEndDate()));
                    }
                }
            }
        }
        System.out.flush();
        System.exit(0);
    }
}
