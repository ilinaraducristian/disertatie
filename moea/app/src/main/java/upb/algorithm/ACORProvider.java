package upb.algorithm;

import org.moeaframework.core.Algorithm;
import org.moeaframework.core.Problem;
import org.moeaframework.core.spi.AlgorithmProvider;
import org.moeaframework.util.TypedProperties;

import java.util.HashMap;
import java.util.Map;
import java.util.Optional;

public class ACORProvider extends AlgorithmProvider {

    private final Map<String, Class<?>> classes = new HashMap<>();

    public ACORProvider() {
        classes.put("1", MoAcoR.class);
        classes.put("2", MOACORSM.class);
        classes.put("3", OMOACOROLD.class);
        classes.put("4", OMOACORV1.class);
        classes.put("5", OMOACORV2.class);
        classes.put("6", OMOACORV3.class);
        classes.put("7", OMOACOR.class);

        classes.put("MoAcoR", MoAcoR.class);
        classes.put("MOACORSM", MOACORSM.class);
        classes.put("OMOACOROLD", OMOACOROLD.class);
        classes.put("OMOACORV1", OMOACORV1.class);
        classes.put("OMOACORV2", OMOACORV2.class);
        classes.put("OMOACORV3", OMOACORV3.class);
        classes.put("OMOACOR", OMOACOR.class);
    }

    @Override
    public Algorithm getAlgorithm(String name, TypedProperties typedProperties, Problem problem) {
        if (!name.toUpperCase().contains("ACOR")) return null;
        String version = typedProperties.getString("version", "");
        return instantiateAlgorithm(version, problem).orElse(instantiateAlgorithm(name, problem).orElse(null));
    }

    private Optional<Algorithm> instantiateAlgorithm(String key, Problem problem) {
        if (!classes.containsKey(key)) {
            return Optional.empty();
        }
        try {
            return Optional.of((Algorithm) classes.get(key).getDeclaredConstructor(Problem.class).newInstance(problem));
        } catch (Exception e) {
            return Optional.empty();
        }
    }

}
