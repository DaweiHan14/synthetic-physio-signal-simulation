import numpy as np
import pandas as pd
from scipy import stats


class UserFeedbackSimulator:
    def __init__(self, seed=42):
        np.random.seed(seed)
        self.params = {
            'users': {
                'ease_of_use': {'baseline': 5.0, 'adaptation': 0.2},
                'comfort': {'baseline': 4.5, 'adaptation': 0.15},
                'reliability': {'baseline': 6.0, 'adaptation': 0.3},
                'satisfaction': {'baseline': 50, 'adaptation': 1.5}
            },
            'controls': {
                'ease_of_use': {'baseline': 7.0, 'adaptation': 0.1},
                'comfort': {'baseline': 6.5, 'adaptation': 0.05},
                'reliability': {'baseline': 8.0, 'adaptation': 0.1},
                'satisfaction': {'baseline': 75, 'adaptation': 0.5}
            }
        }

    def generate_feedback(self, group_size, group_type, weeks=3):
        """Generate simulated user feedback data"""
        data = []
        group_params = self.params[group_type]

        for participant_id in range(1, group_size + 1):
            participant_data = {'participant_id': f"{group_type[:3]}-{participant_id}", 'group': group_type}

            # Weekly feedback data structure
            weekly_scores = {f'week_{week}': {} for week in range(1, weeks + 1)}

            # Generate core metrics
            for metric, params in group_params.items():
                # Baseline + adaptation improvement + individual variability + measurement noise
                trend = [params['baseline'] + params['adaptation'] * week
                         + np.random.normal(0, 0.5)  # individual variability
                         + np.random.normal(0, 0.2)  # measurement noise
                         for week in range(weeks)]

                # Apply logical constraints to keep scores within valid ranges
                if metric == 'satisfaction':
                    trend = np.clip(trend, 0, 100)
                else:
                    trend = np.clip(trend, 1, 10)

                for i, week in enumerate(range(1, weeks + 1)):
                    weekly_scores[f'week_{week}'][metric] = round(trend[i], 1)

            participant_data.update(weekly_scores)
            data.append(participant_data)

        return pd.DataFrame(data)

    def simulate_study(self):
        """Simulate a full study dataset"""
        users_df = self.generate_feedback(group_size=30, group_type='users')
        controls_df = self.generate_feedback(group_size=15, group_type='controls')
        return pd.concat([users_df, controls_df], ignore_index=True)

    def add_qualitative_feedback(self, df):
        """Add simulated qualitative feedback (optional)"""
        # Simplified here; NLG methods could be used for more complex text generation
        positive_comments = [
            "The device is intuitive to use", "Significant improvement in daily activities", "High comfort level",
            "Fast response speed", "Easy to adapt", "Increased independence"
        ]

        negative_comments = [
            "Occasional delays", "Challenging at the beginning", "Some activities remain limited",
            "Needs lighter materials", "Would like easier parameter adjustments", "Occasional noise detected"
        ]

        def generate_comment(row):
            if row['group'] == 'users':
                base_satisfaction = np.mean([row[f'week_{w}']['satisfaction'] for w in [2, 3]])
                if base_satisfaction > 65:
                    return np.random.choice(positive_comments, 2).tolist()
                else:
                    return np.random.choice(negative_comments + positive_comments, 3).tolist()
            return []

        df['qualitative_feedback'] = df.apply(generate_comment, axis=1)
        return df


# Run simulator
if __name__ == "__main__":

    simulator = UserFeedbackSimulator()

    study_data = simulator.simulate_study()

    study_data = simulator.add_qualitative_feedback(study_data)

    # Save results
    study_data.to_csv('simulated_user_feedback.csv', index=False)
    print(study_data.head())

    # Basic validation: check distribution of final satisfaction scores
    for group in ['users', 'controls']:
        group_data = study_data[study_data['group'] == group]
        final_satisfaction = group_data['week_3'].apply(lambda d: d['satisfaction'])

        print(f"{group} final satisfaction scores: μ={final_satisfaction.mean():.1f}, σ={final_satisfaction.std():.1f}")
