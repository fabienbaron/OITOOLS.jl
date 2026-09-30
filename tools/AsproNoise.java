// AsproNoise — compute an ASPRO 2 observation without the GUI and write the OIFITS,
// error bars included. Run under Xvfb: a display must exist, because
// ConfigurationManager skips loading the instrument configuration when headless.
//
//   javac -cp Aspro2-26.09.jar AsproNoise.java
//   xvfb-run -a java -cp .:Aspro2-26.09.jar AsproNoise --list
//   xvfb-run -a java -cp .:Aspro2-26.09.jar AsproNoise --out obs.oifits [options]
//
// PINNED AGAINST Aspro2 26.09 — nothing in CI checks this, so it will break silently:
//
//   jar   : https://www.jmmc.fr/apps/public/Aspro2/Aspro2-26.09.jar   (19 MB)
//           manifest Built-Date 2026/09/10, reports "Aspro2 v26.09" at startup
//   config: CHARA 2026A, bundled in that jar and byte-identical to aspro-conf dd79c85
//           (2026-09-24), src/main/resources/fr/jmmc/aspro/model/CHARA.xml
//   java  : built and run on OpenJDK 11
//
// The signatures this file depends on, which are what a newer release would move:
//
//   Bootstrapper.launchApp(App)
//   ObservationManager: setInterferometerConfigurationName / setInstrumentConfigurationName /
//                       setInstrumentConfigurationStations / setInstrumentMode / setWhen /
//                       setInstrumentSamplingPeriod / setInstrumentAcquisitionTime /
//                       setMinElevation / setAtmosphereQuality / getMainObservation
//   ObservabilityService(ObservationSetting).compute()
//   UVCoverageService(obs, obsData, targetName, uvMax, doUVSupport, useInstrumentBias,
//                     doDataNoise, OIFitsProducerOptions).compute()
//   UVCoverageData.getOiFitsCreator().createOIFits(), OIFitsWriter.writeOIFits(String, OIFitsFile)
//   OIFitsProducerOptions(boolean, boolean, int, UserModelService.MathMode, double)
//
// To move to another release: download that jar, recompile, and run `--list` first — it
// exercises the boot path and the configuration without computing anything.

import fr.jmmc.jmcs.Bootstrapper;
import fr.jmmc.aspro.Aspro2;
import fr.jmmc.aspro.model.ConfigurationManager;
import fr.jmmc.aspro.model.ObservationManager;
import fr.jmmc.aspro.model.oi.ObservationSetting;
import fr.jmmc.aspro.model.oi.Target;
import fr.jmmc.aspro.model.observability.ObservabilityData;
import fr.jmmc.aspro.model.uvcoverage.UVCoverageData;
import fr.jmmc.aspro.service.ObservabilityService;
import fr.jmmc.aspro.service.OIFitsProducerOptions;
import fr.jmmc.aspro.service.UVCoverageService;
import fr.jmmc.aspro.service.UserModelService;
import fr.jmmc.oitools.model.OIFitsFile;
import fr.jmmc.oitools.model.OIFitsWriter;

import java.io.File;
import java.text.SimpleDateFormat;
import java.util.HashMap;
import java.util.Map;

public class AsproNoise {

    static Map<String, String> ARGS = new HashMap<String, String>();

    static String s(String k, String dflt) {
        String v = ARGS.get(k);
        return (v == null) ? dflt : v;
    }

    static double d(String k, double dflt) {
        String v = ARGS.get(k);
        return (v == null) ? dflt : Double.parseDouble(v);
    }

    static Double dObj(String k) {
        String v = ARGS.get(k);
        return (v == null) ? null : Double.valueOf(v);
    }

    static boolean b(String k, boolean dflt) {
        String v = ARGS.get(k);
        return (v == null) ? dflt : Boolean.parseBoolean(v);
    }

    public static void main(String[] argv) throws Exception {
        for (int i = 0; i < argv.length; i++) {
            if (argv[i].startsWith("--")) {
                String key = argv[i].substring(2);
                String val = (i + 1 < argv.length && !argv[i + 1].startsWith("--")) ? argv[++i] : "true";
                ARGS.put(key, val);
            }
        }

        // The app must be booted: the configuration's version check reads ApplicationDescription,
        // and its failure path opens a MODAL dialog that nothing will ever click.
        Bootstrapper.launchApp(new Aspro2(new String[]{}));

        if (ARGS.containsKey("debugNoise")) {
            for (String n : new String[]{"fr.jmmc.aspro.service.NoiseService",
                                         "fr.jmmc.aspro.service.OIFitsCreatorService"}) {
                ch.qos.logback.classic.Logger lg =
                        (ch.qos.logback.classic.Logger) org.slf4j.LoggerFactory.getLogger(n);
                lg.setLevel(ch.qos.logback.classic.Level.DEBUG);
            }
        }

        ConfigurationManager cm = ConfigurationManager.getInstance();
        String interf = s("interferometer", "CHARA");

        if (ARGS.containsKey("list")) {
            System.out.println("\n=== interferometers ===\n  " + cm.getInterferometerNames());
            for (String conf : cm.getInterferometerConfigurationNames(interf)) {
                System.out.println("\n=== configuration: " + conf + " ===");
                for (String ins : cm.getInterferometerInstrumentNames(conf)) {
                    System.out.println("  instrument " + ins);
                    System.out.println("    stations: " + cm.getInstrumentConfigurationNames(conf, ins));
                    System.out.println("    modes   : " + cm.getInstrumentModes(conf, ins));
                }
            }
            System.out.flush();
            System.exit(0);
        }

        String conf     = s("config", cm.getInterferometerConfigurationNames(interf).lastElement());
        String instrum  = s("instrument", "MIRCX-MYSTIC");
        String stations = s("stations", "S1 S2 E1 E2 W1 W2");
        String mode     = s("mode", "Low_H");
        String name     = s("target", "Vega");
        String out      = s("out", "aspro_noise.oifits");

        ObservationManager om = ObservationManager.getInstance();
        om.reset();
        om.setInterferometerConfigurationName(conf);
        om.setInstrumentConfigurationName(instrum);
        om.setInstrumentConfigurationStations(new Object[]{stations});
        om.setInstrumentMode(mode);
        om.setWhen(new SimpleDateFormat("yyyy-MM-dd").parse(s("date", "2026-09-30")));
        // The GUI fills these from the instrument defaults; without them UVCoverageService
        // dereferences a null sampling period.
        om.setInstrumentSamplingPeriod(Double.valueOf(d("sampling", 40.0)));      // minutes between points
        om.setInstrumentAcquisitionTime(Double.valueOf(d("acqTime", 600.0)));     // seconds per point
        om.setMinElevation(d("minElev", 45.0));
        om.setAtmosphereQuality(fr.jmmc.aspro.model.oi.AtmosphereQuality.fromValue(s("atm", "Average")));

        Target t = new Target();
        t.setName(name);
        t.setRA(s("ra", "18:36:56.336"));
        t.setDEC(s("dec", "+38:47:01.28"));
        if (dObj("magV") != null) t.setFLUXV(dObj("magV"));
        if (dObj("magR") != null) t.setFLUXR(dObj("magR"));
        if (dObj("magJ") != null) t.setFLUXJ(dObj("magJ"));
        if (dObj("magH") != null) t.setFLUXH(dObj("magH"));
        if (dObj("magK") != null) t.setFLUXK(dObj("magK"));

        ObservationSetting obs = om.getMainObservation();
        obs.getTargets().add(t);

        System.out.println("\n>>> " + interf + " / " + conf + " / " + instrum
                + " [" + stations + "] mode " + mode + ", target " + name);

        ObservabilityData od = new ObservabilityService(obs).compute();

        OIFitsProducerOptions opt = new OIFitsProducerOptions(
                false, false, (int) d("supersampling", 5),
                UserModelService.MathMode.DEFAULT, d("snrThreshold", 0.0));

        UVCoverageData uvd = new UVCoverageService(obs, od, name,
                d("uvMax", 331.0),
                /* doUVSupport      */ false,
                /* useInstrumentBias*/ b("bias", true),
                /* doDataNoise      */ b("noise", false),
                opt).compute();

        if (uvd == null || uvd.getOiFitsCreator() == null) {
            System.out.println("!!! no OIFits produced — target probably not observable on that date");
            System.exit(2);
        }
        OIFitsFile f = uvd.getOiFitsCreator().createOIFits();
        OIFitsWriter.writeOIFits(new File(out).getAbsolutePath(), f);
        System.out.println(">>> wrote " + out);
        System.out.flush();
        System.exit(0);
    }
}
